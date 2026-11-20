"""Internal helpers of pystackreg.

Nothing in this module is part of the public API; it may change without notice.
"""

from __future__ import annotations

import math
import os
from collections.abc import Callable, Sequence
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from types import TracebackType
from typing import Any

import numpy as np

from pystackreg.util import simple_slice, to_int_dtype

#: Approximate number of array elements processed per chunk by the streaming
#: helpers (16 MiB of float64).
_CHUNK_ELEMENTS = 2**21


def available_cpus() -> int:
    """Number of CPUs usable by the current process (at least 1)."""
    process_cpu_count = getattr(os, "process_cpu_count", None)  # Python >= 3.13
    if process_cpu_count is not None:
        return process_cpu_count() or 1
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0)) or 1
    return os.cpu_count() or 1


def resolve_n_jobs(n_jobs: int | None) -> int:
    """Translate an ``n_jobs`` argument into a number of worker threads.

    Follows the joblib/scikit-learn convention: ``None`` and ``1`` mean
    sequential execution, a positive number is used as is, and a negative
    number ``-k`` means "all CPUs but ``k - 1``" (``-1``: all CPUs).
    """
    if n_jobs is None:
        return 1
    if isinstance(n_jobs, bool) or not isinstance(n_jobs, (int, np.integer)):
        raise TypeError(f"n_jobs must be an integer or None, got {n_jobs!r}")
    n_jobs = int(n_jobs)
    if n_jobs == 0:
        raise ValueError("n_jobs must not be 0; use 1 for sequential execution")
    if n_jobs < 0:
        return max(1, available_cpus() + 1 + n_jobs)
    return n_jobs


def run_tasks(
    func: Callable[[int], None],
    indices: Sequence[int],
    n_jobs: int,
    on_done: Callable[[int], None],
) -> None:
    """Call ``func(i)`` for every index, sequentially or in a thread pool.

    ``func`` must store its own result and must be safe to call concurrently
    for different indices. ``on_done(count)`` is called in the calling thread
    after each completed task with the number of tasks finished so far, so
    progress reporting never runs in a worker thread.

    At most ``2 * n_jobs`` tasks are submitted ahead of time. If a task (or
    ``on_done``) raises, including ``KeyboardInterrupt``, tasks that have not
    started are cancelled, running ones are allowed to finish, and the
    exception is re-raised. Which task's exception is reported first is not
    deterministic when several fail.
    """
    n_workers = min(n_jobs, len(indices))
    if n_workers <= 1:
        for count, index in enumerate(indices, 1):
            func(index)
            on_done(count)
        return

    pending: set[Future[None]] = set()
    remaining = iter(indices)
    finished = 0
    with ThreadPoolExecutor(
        max_workers=n_workers, thread_name_prefix="pystackreg"
    ) as executor:
        try:
            while True:
                for index in remaining:
                    pending.add(executor.submit(func, index))
                    if len(pending) >= 2 * n_workers:
                        break
                if not pending:
                    break
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    future.result()  # re-raises an exception from the worker
                    finished += 1
                    on_done(finished)
        except BaseException:
            for future in pending:
                future.cancel()
            raise


class Progress:
    """Progress reporting through a tqdm bar and/or a user callback.

    :param total: number of iterations
    :param verbose: show a tqdm progress bar
    :param callback: called after every iteration with the keyword arguments
        ``current_iteration`` (1-based count of finished iterations) and
        ``end_iteration`` (the total)
    :param desc: label of the progress bar
    """

    def __init__(
        self,
        total: int,
        verbose: bool,
        callback: Callable[..., Any] | None,
        desc: str,
    ) -> None:
        if callback is not None and not callable(callback):
            raise TypeError("progress_callback must be callable")
        self._total = total
        self._callback = callback
        self._bar = None
        if verbose and total > 0:
            # imported lazily: tqdm.auto selects a notebook widget or a
            # terminal bar and is only needed when a bar is requested
            from tqdm.auto import tqdm

            self._bar = tqdm(total=total, desc=desc)

    def step(self, count: int) -> None:
        """Report that ``count`` iterations have finished."""
        if self._bar is not None:
            self._bar.update(1)
        if self._callback is not None:
            self._callback(current_iteration=count, end_iteration=self._total)

    def __enter__(self) -> Progress:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        if self._bar is not None:
            self._bar.close()
            self._bar = None


class FrameSource:
    """Frame-wise access to a stack, optionally through a moving average.

    The moving average reproduces :func:`pystackreg.util.running_mean`
    (window ``[i + 1 - ceil(N / 2), i + N - ceil(N / 2)]`` with edge
    replication) but computes one frame at a time, so the averaged stack is
    never held in memory.

    :param img: 3-D stack
    :param axis: frame axis (non-negative)
    :param window: moving-average window; 1 disables averaging
    """

    def __init__(self, img: np.ndarray, axis: int, window: int) -> None:
        self._img = img
        self._axis = axis
        self._window = window
        self.n_frames = img.shape[axis]
        if img.dtype.kind == "u":
            self._sum_dtype: Any = np.uint64
        elif img.dtype.kind in "ib":
            self._sum_dtype = np.int64
        else:
            self._sum_dtype = np.float64

    def frame(self, index: int) -> np.ndarray:
        """Return frame ``index`` (a view if no averaging is requested)."""
        if self._window == 1:
            return simple_slice(self._img, index, self._axis)
        first = index + 1 - math.ceil(self._window / 2)
        idx = np.clip(np.arange(first, first + self._window), 0, self.n_frames - 1)
        frames = np.take(self._img, idx, axis=self._axis)
        return frames.sum(axis=self._axis, dtype=self._sum_dtype) / float(self._window)

    def mean_of_first(self, n: int) -> np.ndarray:
        """Mean of the first ``n`` (averaged) frames."""
        if self._window == 1:
            return np.mean(self._img.take(range(n), axis=self._axis), axis=self._axis)
        return np.mean(np.stack([self.frame(i) for i in range(n)]), axis=0)

    def mean_of_all(self) -> np.ndarray:
        """Mean of all (averaged) frames."""
        if self._window == 1:
            return self._img.mean(axis=self._axis)
        total = np.zeros(self.frame(0).shape, dtype=np.float64)
        for i in range(self.n_frames):
            total += self.frame(i)
        return total / self.n_frames


def detect_time_axis(img: np.ndarray) -> int:
    """Return the axis of a 3-D stack with the lowest average variance.

    Equivalent to ``argmin([np.var(img, axis=i).mean() for i in range(3)])``,
    but evaluated in chunks along the first axis: ``np.var`` on the whole
    stack allocates a float64 copy of it, which for integer stacks is several
    times the size of the data.
    """
    n0 = img.shape[0]
    if img.size == 0:
        return 0
    plane = img.shape[1] * img.shape[2]
    step = max(1, _CHUNK_ELEMENTS // max(plane, 1))

    # Variance along axis 0 from running sums of values shifted by the first
    # plane (the shift keeps the sums well conditioned for data with a large
    # offset). Variances along axes 1 and 2 are local to a chunk.
    first = img[0].astype(np.float64)
    sum0 = np.zeros_like(first)
    sumsq0 = np.zeros_like(first)
    var1_total = 0.0
    var2_total = 0.0
    for start in range(0, n0, step):
        chunk = img[start : start + step].astype(np.float64)
        var1_total += float(np.var(chunk, axis=1).sum())
        var2_total += float(np.var(chunk, axis=2).sum())
        chunk -= first
        sum0 += chunk.sum(axis=0)
        sumsq0 += np.square(chunk, out=chunk).sum(axis=0)
    var0 = float((sumsq0 / n0 - np.square(sum0 / n0)).mean())
    var1 = var1_total / (n0 * img.shape[2])
    var2 = var2_total / (n0 * img.shape[1])
    return int(np.argmin([var0, var1, var2]))


def check_output_array(out: Any, shape: tuple[int, ...]) -> np.ndarray:
    """Validate a user-supplied output array for :func:`store_frame`."""
    if not isinstance(out, np.ndarray):
        raise TypeError("out must be a numpy array (np.memmap is accepted)")
    if out.shape != shape:
        raise ValueError(f"out has shape {out.shape}, but the stack has shape {shape}")
    if out.dtype.kind not in "fiu":
        raise TypeError(
            f"out must have a floating-point or integer dtype, got {out.dtype}"
        )
    if not out.flags.writeable:
        raise ValueError("out is read-only")
    return out


def store_frame(out: np.ndarray, index: int, axis: int, frame: np.ndarray) -> None:
    """Write a transformed (float64) frame into ``out`` along ``axis``.

    Integer outputs receive values rounded and clipped to the range of the
    dtype (see :func:`pystackreg.util.to_int_dtype`); a plain cast would
    truncate and wrap around.
    """
    target: Any = [slice(None)] * out.ndim
    target[axis] = index
    if out.dtype.kind in "iu":
        frame = to_int_dtype(frame, out.dtype.type)
    out[tuple(target)] = frame
