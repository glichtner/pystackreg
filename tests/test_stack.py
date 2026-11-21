"""Stack registration and transformation: references, axes, moving average,
memory-related options and state handling."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from pystackreg import StackReg
from pystackreg import _internal as internal
from pystackreg.util import running_mean, to_uint16

from helpers import TRANSFORMATIONS, drifting_stack


@pytest.fixture(scope="module")
def stack():
    data, _ = drifting_stack(n_frames=6, shape=(48, 60), kind="rigid_body")
    return np.round(data).astype(np.uint16)


@pytest.fixture(scope="module")
def tmats(stack):
    return StackReg(StackReg.RIGID_BODY).register_stack(stack, reference="previous")


# -- references ---------------------------------------------------------------


def test_previous_chains_pairwise_registrations(stack, tmats):
    sr = StackReg(StackReg.RIGID_BODY)
    expected = [np.identity(3)]
    for i in range(1, len(stack)):
        expected.append(sr.register(stack[i - 1], stack[i]) @ expected[-1])
    np.testing.assert_array_equal(tmats, np.array(expected))


def test_first_registers_every_frame_to_frame_0(stack):
    sr = StackReg(StackReg.RIGID_BODY)
    result = sr.register_stack(stack, reference="first")
    np.testing.assert_array_equal(result[0], np.identity(3))
    for i in range(1, len(stack)):
        np.testing.assert_array_equal(result[i], sr.register(stack[0], stack[i]))


def test_first_with_n_frames_uses_their_mean(stack):
    sr = StackReg(StackReg.RIGID_BODY)
    result = sr.register_stack(stack, reference="first", n_frames=3)
    ref = stack[:3].mean(axis=0)
    np.testing.assert_array_equal(result[0], np.identity(3))
    np.testing.assert_array_equal(result[4], sr.register(ref, stack[4]))


def test_mean_registers_every_frame_to_the_mean_image(stack):
    sr = StackReg(StackReg.RIGID_BODY)
    result = sr.register_stack(stack, reference="mean")
    ref = stack.mean(axis=0)
    for i in range(len(stack)):
        np.testing.assert_array_equal(result[i], sr.register(ref, stack[i]))


def test_bilinear_stack_with_fixed_reference(stack):
    result = StackReg(StackReg.BILINEAR).register_stack(stack, reference="first")
    assert result.shape == (len(stack), 4, 4)
    np.testing.assert_array_equal(result[0], np.identity(4))


# -- frame axis ---------------------------------------------------------------


@pytest.mark.parametrize("axis", [1, 2, -1, -2, -3])
@pytest.mark.parametrize("reference", ["previous", "first", "mean"])
def test_frame_axis(stack, axis, reference):
    sr = StackReg(StackReg.RIGID_BODY)
    expected_tmats = sr.register_stack(stack, reference=reference)
    expected = sr.transform_stack(stack)

    moved = np.moveaxis(stack, 0, axis)
    sr = StackReg(StackReg.RIGID_BODY)
    result = sr.register_transform_stack(moved, reference=reference, axis=axis)
    assert result.shape == moved.shape
    np.testing.assert_array_equal(sr._tmats, expected_tmats)
    np.testing.assert_array_equal(np.moveaxis(result, axis, 0), expected)


def test_axis_warning(stack):
    sr = StackReg(StackReg.TRANSLATION)
    with pytest.warns(UserWarning, match="possible time axis"):
        sr.register_stack(stack, axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        sr.register_stack(stack, axis=1, suppress_axis_warning=True)
        sr.register_transform_stack(stack, axis=1, suppress_axis_warning=True)
        sr.register_stack(stack, axis=0)


@pytest.mark.parametrize("shape", [(7, 9, 11), (3, 20, 4), (12, 5, 5)])
@pytest.mark.parametrize("dtype", ["uint16", "float64"])
def test_time_axis_detection_matches_definition(shape, dtype, monkeypatch):
    """The chunked evaluation equals argmin of the mean variance per axis."""
    rng = np.random.default_rng(0)
    for axis in range(3):
        scale = np.ones(3)
        scale[axis] = 0.05  # little variation along this axis
        grids = np.meshgrid(*[np.arange(n) for n in shape], indexing="ij")
        data = 500 + sum(
            200 * s * np.sin(g * 1.3 + i)
            for i, (s, g) in enumerate(zip(scale, grids, strict=True))
        )
        data = (data + rng.normal(0, 1, shape)).astype(dtype)
        expected = int(np.argmin([np.var(data, axis=i).mean() for i in range(3)]))
        assert expected == axis
        assert StackReg._detect_time_axis(data) == expected
        monkeypatch.setattr(internal, "_CHUNK_ELEMENTS", 50)  # force many chunks
        assert StackReg._detect_time_axis(data) == expected
        monkeypatch.undo()


# -- moving average -----------------------------------------------------------


@pytest.mark.parametrize("window", [2, 3, 4])
@pytest.mark.parametrize("reference", ["first", "mean"])
def test_moving_average_equals_running_mean(stack, window, reference):
    """The frame-wise moving average is the one of ``util.running_mean``."""
    sr = StackReg(StackReg.RIGID_BODY)
    result = sr.register_stack(stack, reference=reference, moving_average=window)
    averaged = running_mean(stack, window, axis=0)
    ref = averaged[0] if reference == "first" else averaged.mean(axis=0)
    for i in range(len(stack)):
        np.testing.assert_allclose(
            result[i], sr.register(ref, averaged[i]), rtol=0, atol=1e-9
        )


@pytest.mark.parametrize("window", [0, 1])
def test_moving_average_of_zero_or_one_is_no_average(stack, window):
    """0 switches the moving average off, as 1 does (0.2.x accepted it)."""
    sr = StackReg(StackReg.RIGID_BODY)
    expected = sr.register_stack(stack, reference="first")
    result = sr.register_stack(stack, reference="first", moving_average=window)
    np.testing.assert_array_equal(result, expected)


def test_moving_average_with_previous_keeps_first_frame_fixed(stack):
    """Regression test: 0.2.8 registered frame 0 against the *last* frame."""
    sr = StackReg(StackReg.RIGID_BODY)
    result = sr.register_stack(stack, reference="previous", moving_average=3)
    np.testing.assert_array_equal(result[0], np.identity(3))
    averaged = running_mean(stack, 3, axis=0)
    np.testing.assert_allclose(
        result[1], sr.register(averaged[0], averaged[1]), rtol=0, atol=1e-9
    )


@pytest.mark.parametrize("axis", [1, 2])
def test_moving_average_along_other_axis(stack, axis):
    """Regression test: 0.2.8 produced garbage for axis != 0."""
    expected = StackReg(StackReg.RIGID_BODY).register_stack(
        stack, reference="first", moving_average=3
    )
    result = StackReg(StackReg.RIGID_BODY).register_stack(
        np.moveaxis(stack, 0, axis), reference="first", moving_average=3, axis=axis
    )
    np.testing.assert_array_equal(result, expected)


def test_moving_average_float_stack(stack):
    data = stack.astype(np.float32)
    sr = StackReg(StackReg.RIGID_BODY)
    result = sr.register_stack(data, reference="first", moving_average=3)
    expected = sr.register_stack(stack, reference="first", moving_average=3)
    np.testing.assert_allclose(result, expected, rtol=0, atol=1e-6)


@pytest.mark.parametrize("axis", [0, 1, 2])
@pytest.mark.parametrize("window", [1, 2, 3, 5])
def test_running_mean(axis, window):
    x = np.arange(4 * 5 * 6, dtype=float).reshape(4, 5, 6) ** 1.5
    result = running_mean(x, window, axis=axis)
    assert result.shape == x.shape
    n = x.shape[axis]
    for i in range(n):
        first = i + 1 - int(np.ceil(window / 2))
        idx = np.clip(np.arange(first, first + window), 0, n - 1)
        expected = np.take(x, idx, axis=axis).mean(axis=axis)
        np.testing.assert_allclose(np.take(result, i, axis=axis), expected)


# -- special stacks -----------------------------------------------------------


@pytest.mark.parametrize("reference", ["previous", "first", "mean"])
def test_single_frame_stack(stack, reference):
    sr = StackReg(StackReg.RIGID_BODY)
    single = stack[:1]
    result = sr.register_stack(single, reference=reference, suppress_axis_warning=True)
    assert result.shape == (1, 3, 3)
    np.testing.assert_allclose(result[0], np.identity(3), atol=1e-6)
    out = sr.transform_stack(single)
    np.testing.assert_allclose(out, single, rtol=1e-5, atol=0.05)


def test_two_frame_stack(stack):
    sr = StackReg(StackReg.RIGID_BODY)
    result = sr.register_stack(stack[:2], reference="previous")
    np.testing.assert_array_equal(result[1], sr.register(stack[0], stack[1]))


# -- state --------------------------------------------------------------------


def test_state_after_register_stack(stack):
    sr = StackReg(StackReg.RIGID_BODY)
    assert not sr.is_registered()
    result = sr.register_stack(stack, reference="first")
    assert sr.is_registered()
    # as before 0.3.0, the single-image state reflects the last frame
    np.testing.assert_array_equal(sr.get_matrix(), result[-1])
    refpts, movpts = sr.get_points()
    assert refpts.shape == movpts.shape == (3, 2)
    np.testing.assert_array_equal(
        sr.transform(stack[-1]), sr.transform_stack(stack)[-1]
    )


def test_register_transform_stack_equals_two_steps(stack):
    sr = StackReg(StackReg.AFFINE)
    one = sr.register_transform_stack(stack, reference="first")
    sr2 = StackReg(StackReg.AFFINE)
    sr2.register_stack(stack, reference="first")
    np.testing.assert_array_equal(one, sr2.transform_stack(stack))


def test_saved_matrices_can_be_reused(stack, tmats, tmp_path):
    expected = StackReg(StackReg.RIGID_BODY).transform_stack(stack, tmats=tmats)
    np.save(tmp_path / "tmats.npy", tmats)
    loaded = np.load(tmp_path / "tmats.npy")
    # a fresh object that never registered anything
    result = StackReg(StackReg.RIGID_BODY).transform_stack(stack, tmats=loaded)
    np.testing.assert_array_equal(result, expected)
    # matrices can be applied to another stack of the same length (e.g.
    # another colour channel)
    other = stack[:, ::-1, :].copy()
    out = StackReg(StackReg.RIGID_BODY).transform_stack(other, tmats=loaded)
    assert out.shape == other.shape


def test_transform_stack_default_output(stack, tmats):
    out = StackReg(StackReg.RIGID_BODY).transform_stack(stack, tmats=tmats)
    assert out.dtype == np.float64
    assert out.flags.c_contiguous
    np.testing.assert_array_equal(stack[0], np.round(out[0]))


# -- out= -----------------------------------------------------------------------


def test_out_float32(stack, tmats):
    sr = StackReg(StackReg.RIGID_BODY)
    expected = sr.transform_stack(stack, tmats=tmats)
    out = np.empty(stack.shape, dtype=np.float32)
    result = sr.transform_stack(stack, tmats=tmats, out=out)
    assert result is out
    # the default output already holds single-precision values
    np.testing.assert_array_equal(out, expected)


def test_out_integer_rounds_and_clips(stack, tmats):
    sr = StackReg(StackReg.RIGID_BODY)
    expected = to_uint16(sr.transform_stack(stack, tmats=tmats))
    out = np.empty(stack.shape, dtype=np.uint16)
    sr.transform_stack(stack, tmats=tmats, out=out)
    np.testing.assert_array_equal(out, expected)

    small = np.empty(stack.shape, dtype=np.uint8)
    sr.transform_stack(stack, tmats=tmats, out=small)
    assert small.max() == 255  # clipped, not wrapped around


def test_out_in_place(stack, tmats):
    sr = StackReg(StackReg.RIGID_BODY)
    expected = to_uint16(sr.transform_stack(stack, tmats=tmats))
    work = stack.copy()
    result = sr.transform_stack(work, tmats=tmats, out=work)
    assert result is work
    np.testing.assert_array_equal(work, expected)


def test_out_memmap(stack, tmats, tmp_path):
    sr = StackReg(StackReg.RIGID_BODY)
    expected = sr.transform_stack(stack, tmats=tmats)
    out = np.lib.format.open_memmap(
        tmp_path / "registered.npy", mode="w+", dtype=np.float32, shape=stack.shape
    )
    sr.transform_stack(stack, tmats=tmats, out=out)
    out.flush()
    del out
    np.testing.assert_array_equal(np.load(tmp_path / "registered.npy"), expected)


@pytest.mark.parametrize("axis", [1, 2])
def test_out_with_other_axis(stack, tmats, axis):
    sr = StackReg(StackReg.RIGID_BODY)
    moved = np.moveaxis(stack, 0, axis)
    out = np.empty(moved.shape, dtype=np.float32)
    sr.transform_stack(moved, axis=axis, tmats=tmats, out=out)
    np.testing.assert_array_equal(
        np.moveaxis(out, axis, 0), sr.transform_stack(stack, tmats=tmats)
    )


def test_register_transform_stack_out(stack):
    sr = StackReg(StackReg.RIGID_BODY)
    expected = to_uint16(sr.register_transform_stack(stack, reference="first"))
    out = np.empty_like(stack)
    result = sr.register_transform_stack(stack, reference="first", out=out)
    assert result is out
    np.testing.assert_array_equal(out, expected)


# -- input variants -------------------------------------------------------------


@pytest.mark.parametrize("name", list(TRANSFORMATIONS))
def test_integer_and_float_stacks_agree(stack, name):
    reference = "first" if name == "bilinear" else "previous"
    a = StackReg(TRANSFORMATIONS[name]).register_stack(stack, reference=reference)
    b = StackReg(TRANSFORMATIONS[name]).register_stack(
        stack.astype(np.float64), reference=reference
    )
    np.testing.assert_array_equal(a, b)


def test_stack_layout_does_not_matter(stack, tmats):
    sr = StackReg(StackReg.RIGID_BODY)
    expected = sr.transform_stack(stack, tmats=tmats)

    readonly = stack.copy()
    readonly.setflags(write=False)
    strided = np.repeat(stack, 2, axis=2)[:, :, ::2]
    assert not strided.flags.c_contiguous
    for variant in (np.asfortranarray(stack), readonly, strided, stack.tolist()):
        result = StackReg(StackReg.RIGID_BODY).register_stack(variant)
        np.testing.assert_array_equal(result, tmats)
        np.testing.assert_array_equal(
            sr.transform_stack(variant, tmats=tmats), expected
        )
    np.testing.assert_array_equal(readonly, stack)


def test_memmap_input(stack, tmats, tmp_path):
    np.save(tmp_path / "stack.npy", stack)
    mapped = np.load(tmp_path / "stack.npy", mmap_mode="r")
    result = StackReg(StackReg.RIGID_BODY).register_stack(mapped)
    np.testing.assert_array_equal(result, tmats)
