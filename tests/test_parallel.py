"""Parallel execution, thread safety and progress reporting."""

from __future__ import annotations

import threading
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

from pystackreg import StackReg, turboreg
from pystackreg import _internal as internal

from helpers import TRANSFORMATIONS, Scene, drifting_stack, known_transformations


@pytest.fixture(scope="module")
def stack():
    data, _ = drifting_stack(n_frames=10, shape=(48, 60), kind="rigid_body")
    return np.round(data).astype(np.uint16)


# -- equivalence with sequential execution -------------------------------------


@pytest.mark.parametrize("n_jobs", [2, 4, -1])
@pytest.mark.parametrize("reference", ["previous", "first", "mean"])
@pytest.mark.parametrize("name", list(TRANSFORMATIONS))
def test_parallel_registration_equals_sequential(stack, name, reference, n_jobs):
    if name == "bilinear" and reference == "previous":
        pytest.skip("bilinear transformations cannot be chained")
    sequential = StackReg(TRANSFORMATIONS[name])
    expected = sequential.register_stack(stack, reference=reference)
    parallel = StackReg(TRANSFORMATIONS[name])
    result = parallel.register_stack(stack, reference=reference, n_jobs=n_jobs)
    # bit-identical: every frame is an independent computation
    np.testing.assert_array_equal(result, expected)
    np.testing.assert_array_equal(parallel.get_matrix(), sequential.get_matrix())
    for a, b in zip(parallel.get_points(), sequential.get_points(), strict=True):
        np.testing.assert_array_equal(a, b)


@pytest.mark.parametrize("n_jobs", [2, 4, -1])
def test_parallel_transformation_equals_sequential(stack, n_jobs):
    sr = StackReg(StackReg.AFFINE)
    sr.register_stack(stack, reference="first")
    expected = sr.transform_stack(stack)
    np.testing.assert_array_equal(sr.transform_stack(stack, n_jobs=n_jobs), expected)
    out = np.empty(stack.shape, dtype=np.uint16)
    sr.transform_stack(stack, n_jobs=n_jobs, out=out)
    np.testing.assert_array_equal(out, np.round(expected.clip(0, 65535)))


def test_parallel_with_moving_average_and_axis(stack):
    moved = np.moveaxis(stack, 0, 2)
    kwargs = {"reference": "previous", "moving_average": 3, "axis": 2}
    expected = StackReg(StackReg.RIGID_BODY).register_transform_stack(moved, **kwargs)
    result = StackReg(StackReg.RIGID_BODY).register_transform_stack(
        moved, n_jobs=3, **kwargs
    )
    np.testing.assert_array_equal(result, expected)


def test_repeated_parallel_runs_are_deterministic(stack):
    runs = [
        StackReg(StackReg.RIGID_BODY).register_stack(stack, n_jobs=4) for _ in range(5)
    ]
    for run in runs[1:]:
        np.testing.assert_array_equal(run, runs[0])


def test_more_workers_than_frames(stack):
    expected = StackReg(StackReg.RIGID_BODY).register_stack(stack[:3])
    result = StackReg(StackReg.RIGID_BODY).register_stack(stack[:3], n_jobs=32)
    np.testing.assert_array_equal(result, expected)


def test_n_jobs_resolution(monkeypatch):
    monkeypatch.setattr(internal, "available_cpus", lambda: 8)
    assert internal.resolve_n_jobs(None) == 1
    assert internal.resolve_n_jobs(1) == 1
    assert internal.resolve_n_jobs(np.int64(3)) == 3
    assert internal.resolve_n_jobs(-1) == 8
    assert internal.resolve_n_jobs(-2) == 7
    assert internal.resolve_n_jobs(-100) == 1
    assert internal.available_cpus() >= 1


# -- thread safety of the native code -----------------------------------------


def test_concurrent_calls_from_independent_threads():
    """Many threads use the extension at the same time, each with its own
    StackReg object and different data; every result must equal the result
    of the same call made alone."""
    jobs = []
    for seed, name in enumerate(list(TRANSFORMATIONS) * 4):
        shape = (40 + 4 * (seed % 5), 56 + 2 * (seed % 3))
        scene = Scene(shape, seed=seed + 1)
        jobs.append(
            (name, scene.render(), scene.render(known_transformations(shape)[name]))
        )

    def work(job):
        name, ref, mov = job
        sr = StackReg(TRANSFORMATIONS[name])
        tmat = sr.register(ref, mov)
        return tmat, sr.transform(mov)

    expected = [work(job) for job in jobs]
    with ThreadPoolExecutor(max_workers=8) as executor:
        for _ in range(3):
            results = list(executor.map(work, jobs))
            for (tmat, out), (exp_tmat, exp_out) in zip(results, expected, strict=True):
                np.testing.assert_array_equal(tmat, exp_tmat)
                np.testing.assert_array_equal(out, exp_out)


def test_native_code_releases_the_gil():
    """While one thread is inside the extension, another thread must be able
    to execute Python code.

    A background thread appends a timestamp about every half millisecond, but
    only while the main thread is inside a *single* native call. If the
    extension held the GIL, the background thread could not run during the
    call at all (at most once, right after the call returns); with the GIL
    released it runs hundreds of times.
    """
    scene = Scene((640, 640), n_blobs=200, seed=9)
    ref = scene.render()
    mov = scene.render(known_transformations(scene.shape)["affine"])

    ticks = []
    inside = threading.Event()
    stop = threading.Event()

    def ticker():
        while not stop.is_set():
            if inside.is_set():
                ticks.append(time.perf_counter())
            time.sleep(0.0005)

    def ticks_during(call):
        del ticks[:]
        inside.set()
        start = time.perf_counter()
        call()
        elapsed = time.perf_counter() - start
        inside.clear()
        return len(ticks), elapsed

    thread = threading.Thread(target=ticker)
    thread.start()
    try:
        time.sleep(0.01)
        matrix = np.array([[0.5, 1.0, 0.01], [-0.5, -0.01, 1.0]])
        results = {
            "_register": ticks_during(
                lambda: turboreg._register(
                    ref[:-1, :-1], mov[:-1, :-1], StackReg.AFFINE
                )
            ),
            "_transform": ticks_during(lambda: turboreg._transform(mov, matrix)),
        }
    finally:
        stop.set()
        thread.join()
    for name, (count, elapsed) in results.items():
        assert elapsed > 0.005, f"{name}: workload too small to be meaningful"
        assert count >= 5, f"{name}: only {count} ticks in {elapsed * 1e3:.0f} ms"


# -- errors and cancellation ---------------------------------------------------


@pytest.mark.parametrize("n_jobs", [1, 4])
def test_worker_exception_propagates(stack, n_jobs):
    bad = stack.astype(np.float64)
    bad[5, 2, 2] = np.nan
    sr = StackReg(StackReg.RIGID_BODY)
    with pytest.raises(ValueError, match="NaN or infinite"):
        sr.register_stack(bad, reference="first", n_jobs=n_jobs)
    # the object stays usable
    np.testing.assert_array_equal(
        sr.register_stack(stack, n_jobs=n_jobs),
        StackReg(StackReg.RIGID_BODY).register_stack(stack),
    )
    assert not [t for t in threading.enumerate() if t.name.startswith("pystackreg")]


def test_exception_cancels_outstanding_tasks():
    calls = []

    def task(i):
        calls.append(i)
        if i == 0:
            raise RuntimeError("boom")
        time.sleep(0.02)

    with pytest.raises(RuntimeError, match="boom"):
        internal.run_tasks(task, range(200), 2, lambda count: None)
    # only the tasks queued ahead (2 * n_jobs) and the few started before the
    # failure was noticed have run
    assert len(calls) < 20


@pytest.mark.parametrize("n_jobs", [1, 3])
def test_exception_in_callback_stops_processing(stack, n_jobs):
    seen = []

    def callback(current_iteration, end_iteration):
        seen.append(current_iteration)
        if current_iteration == 3:
            raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        StackReg(StackReg.RIGID_BODY).register_stack(
            stack, n_jobs=n_jobs, progress_callback=callback
        )
    assert seen == [1, 2, 3]
    assert not [t for t in threading.enumerate() if t.name.startswith("pystackreg")]


# -- progress reporting ----------------------------------------------------------


@pytest.mark.parametrize("n_jobs", [1, 4])
@pytest.mark.parametrize(
    ("reference", "skipped"), [("previous", 1), ("first", 1), ("mean", 0)]
)
def test_progress_callback(stack, reference, skipped, n_jobs):
    calls = []
    main_thread = threading.get_ident()

    def callback(current_iteration, end_iteration):
        calls.append((current_iteration, end_iteration, threading.get_ident()))

    StackReg(StackReg.RIGID_BODY).register_stack(
        stack, reference=reference, n_jobs=n_jobs, progress_callback=callback
    )
    total = len(stack) - skipped
    assert [c[0] for c in calls] == list(range(1, total + 1))
    assert {c[1] for c in calls} == {total}
    assert {c[2] for c in calls} == {main_thread}


@pytest.mark.parametrize("n_jobs", [1, 4])
def test_transform_stack_progress_callback(stack, n_jobs):
    sr = StackReg(StackReg.RIGID_BODY)
    sr.register_stack(stack)
    calls = []
    sr.transform_stack(
        stack, n_jobs=n_jobs, progress_callback=lambda **kw: calls.append(kw)
    )
    assert [c["current_iteration"] for c in calls] == list(range(1, len(stack) + 1))
    assert {c["end_iteration"] for c in calls} == {len(stack)}


def test_register_transform_stack_reports_registration_only(stack):
    """Unchanged from 0.2.x: the callback counts registration iterations."""
    calls = []
    StackReg(StackReg.RIGID_BODY).register_transform_stack(
        stack, progress_callback=lambda **kw: calls.append(kw)
    )
    assert len(calls) == len(stack) - 1


def test_silent_by_default(stack, capsys):
    StackReg(StackReg.RIGID_BODY).register_transform_stack(stack, n_jobs=2)
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


@pytest.mark.parametrize("n_jobs", [1, 4])
def test_verbose_shows_progress_bars(stack, capsys, n_jobs):
    tqdm = pytest.importorskip("tqdm")
    StackReg(StackReg.RIGID_BODY).register_transform_stack(
        stack, verbose=True, n_jobs=n_jobs
    )
    err = capsys.readouterr().err
    assert "Registering" in err
    assert "Transforming" in err
    assert f"{len(stack) - 1}/{len(stack) - 1}" in err  # registration completed
    assert f"{len(stack)}/{len(stack)}" in err  # transformation completed
    assert not getattr(tqdm.tqdm, "_instances", ())  # all bars were closed


def test_progress_bar_is_closed_on_error(stack, capsys):
    tqdm = pytest.importorskip("tqdm")
    bad = stack.astype(np.float64)
    bad[4, 2, 2] = np.inf
    with pytest.raises(ValueError, match="NaN or infinite"):
        StackReg(StackReg.RIGID_BODY).register_stack(bad, verbose=True)
    capsys.readouterr()
    assert not getattr(tqdm.tqdm, "_instances", ())
