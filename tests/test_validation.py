"""Defined behaviour for invalid input."""

from __future__ import annotations

import numpy as np
import pytest

from pystackreg import StackReg, turboreg

from helpers import TRANSFORMATIONS, Scene


@pytest.fixture(scope="module")
def pair():
    scene = Scene((40, 50))
    ref = scene.render()
    return ref, np.roll(ref, 1, axis=1)


@pytest.mark.parametrize("value", [0, 1, 5, 7, -1, "affine", None, 2.5])
def test_invalid_transformation(value):
    with pytest.raises(ValueError, match="Invalid transformation"):
        StackReg(value)


def test_errors_are_still_instances_of_exception():
    """Code written against 0.2.x catches ``Exception``."""
    with pytest.raises(Exception):  # noqa: B017
        StackReg(99)


@pytest.mark.parametrize("name", list(TRANSFORMATIONS))
@pytest.mark.parametrize("shape", [(1, 1), (1, 10), (2, 2), (2, 10), (10, 2)])
def test_register_rejects_tiny_images(name, shape):
    img = np.ones(shape)
    with pytest.raises(ValueError, match="at least 3 x 3"):
        StackReg(TRANSFORMATIONS[name]).register(img, img)


@pytest.mark.parametrize("name", list(TRANSFORMATIONS))
@pytest.mark.parametrize(
    "shape", [(3, 3), (3, 9), (5, 4), (13, 12), (24, 24), (25, 47)]
)
def test_register_accepts_small_images(name, shape):
    """Smallest sizes and the sizes around the first pyramid level (24 px)."""
    y, x = np.mgrid[0 : shape[0], 0 : shape[1]].astype(float)
    img = 100 + 50 * np.sin(0.4 * x) * np.cos(0.3 * y) + x + 0.5 * y
    tmat = StackReg(TRANSFORMATIONS[name]).register(img, img)
    assert np.isfinite(tmat).all()


def test_register_rejects_wrong_dimensions(pair):
    ref, mov = pair
    sr = StackReg(StackReg.RIGID_BODY)
    with pytest.raises(ValueError, match="two-dimensional"):
        sr.register(ref[np.newaxis], mov[np.newaxis])
    with pytest.raises(ValueError, match="two-dimensional"):
        sr.register(ref[0], mov[0])
    with pytest.raises(ValueError, match="same shape"):
        sr.register(ref, mov[:, :-1])


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
@pytest.mark.parametrize("which", ["ref", "mov"])
def test_register_rejects_non_finite_values(pair, bad, which):
    ref, mov = (a.copy() for a in pair)
    (ref if which == "ref" else mov)[7, 9] = bad
    with pytest.raises(ValueError, match=f"{which} contains NaN or infinite"):
        StackReg(StackReg.AFFINE).register(ref, mov)


def test_register_rejects_complex_data(pair):
    ref, mov = pair
    with pytest.raises(TypeError):
        StackReg(StackReg.AFFINE).register(ref.astype(complex), mov.astype(complex))


@pytest.mark.parametrize(
    "name", ["translation", "rigid_body", "scaled_rotation", "affine"]
)
def test_constant_images_give_identity(name):
    """Without any structure the optimiser cannot move: the result is the
    identity. (The image gradient and hence the Hessian are exactly zero; the
    native code handles the resulting non-finite updates like the Java
    original.)"""
    flat = np.full((40, 50), 7.0)
    tmat = StackReg(TRANSFORMATIONS[name]).register(flat, flat)
    np.testing.assert_allclose(tmat, np.identity(3), atol=1e-12)


def test_constant_images_bilinear_result_is_undefined_but_finite():
    """The bilinear optimiser differentiates the interpolated image, whose
    gradient is zero only up to rounding noise. For structureless images it
    therefore follows that noise: the result is meaningless, but the call
    must neither fail nor return non-finite values. All-zero images have an
    exactly vanishing gradient and give the identity."""
    flat = np.full((40, 50), 7.0)
    assert np.isfinite(StackReg(StackReg.BILINEAR).register(flat, flat)).all()
    zeros = np.zeros((40, 50))
    np.testing.assert_allclose(
        StackReg(StackReg.BILINEAR).register(zeros, zeros), np.identity(4), atol=1e-12
    )


def test_methods_requiring_a_registration(pair):
    ref, _ = pair
    sr = StackReg(StackReg.AFFINE)
    with pytest.raises(RuntimeError, match="register"):
        sr.get_matrix()
    with pytest.raises(RuntimeError, match="register"):
        sr.transform(ref)
    with pytest.raises(RuntimeError, match="register"):
        sr.transform_points([[1.0, 2.0]])
    with pytest.raises(RuntimeError, match="register"):
        sr.transform_stack(np.stack([ref, ref]))
    assert sr.get_points() == (None, None)


@pytest.mark.parametrize("shape", [(2, 3), (4, 4), (3,), (3, 3, 3)])
def test_set_matrix_rejects_wrong_shape(shape):
    with pytest.raises(ValueError, match="Invalid shape"):
        StackReg(StackReg.AFFINE).set_matrix(np.zeros(shape))


def test_bilinear_needs_4x4_matrix(pair):
    ref, _ = pair
    with pytest.raises(ValueError, match="Invalid shape"):
        StackReg(StackReg.BILINEAR).transform(ref, np.identity(3))
    with pytest.raises(ValueError, match="Invalid shape"):
        StackReg(StackReg.AFFINE).transform(ref, np.identity(4))


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, 1e300, 3e9])
@pytest.mark.parametrize("name", ["translation", "affine", "bilinear"])
def test_transform_rejects_degenerate_matrices(pair, name, bad):
    ref, _ = pair
    tmat = np.identity(4 if name == "bilinear" else 3)
    tmat[0, -1] = bad
    with pytest.raises(ValueError, match="finite"):
        StackReg(TRANSFORMATIONS[name]).transform(ref, tmat)


def test_transform_rejects_wrong_image_dimensions(pair):
    ref, _ = pair
    sr = StackReg(StackReg.AFFINE)
    with pytest.raises(ValueError, match="two dimensional"):
        sr.transform(ref[np.newaxis], np.identity(3))
    with pytest.raises(ValueError, match="at least 1 x 1"):
        sr.transform(np.zeros((0, 5)), np.identity(3))


def test_native_functions_validate_their_arguments(pair):
    """The extension module must not crash on arguments that bypass the
    Python layer."""
    ref, mov = pair
    with pytest.raises(ValueError, match="Invalid transformation"):
        turboreg._register(ref, mov, 5)
    with pytest.raises(ValueError, match="at least 2 x 2"):
        turboreg._register(ref[:1], mov[:1], 3)
    with pytest.raises(ValueError, match="same shape"):
        turboreg._register(ref, mov[:-1], 3)
    with pytest.raises(ValueError, match="two dimensional"):
        turboreg._register(ref[0], mov[0], 3)
    with pytest.raises(ValueError, match=r"shape \(2,1\), \(2,3\) or \(2,4\)"):
        turboreg._transform(ref, np.zeros((3, 3)))
    with pytest.raises((TypeError, ValueError)):
        turboreg._register("not an image", mov, 3)
    with pytest.raises(TypeError):
        turboreg._register(ref)


def test_register_returns_new_references(pair):
    """Guards against the reference leak fixed in 0.3.0: the returned arrays
    must be owned by the caller alone."""
    import sys

    ref, mov = pair
    result = turboreg._register(ref[:-1, :-1], mov[:-1, :-1], StackReg.AFFINE)
    assert isinstance(result, tuple)
    assert len(result) == 3
    for item in result:
        # one reference from the tuple, one from the loop variable, one from
        # the argument of getrefcount
        assert sys.getrefcount(item) == 3


# -- stack functions ---------------------------------------------------------


@pytest.fixture(scope="module")
def stack():
    scene = Scene((32, 40))
    ref = scene.render()
    return np.stack([np.roll(ref, i, axis=1) for i in range(4)])


def test_stack_must_be_three_dimensional(stack):
    sr = StackReg(StackReg.TRANSLATION)
    for bad in (stack[0], stack[np.newaxis]):
        with pytest.raises(ValueError, match="three dimensions"):
            sr.register_stack(bad)
        with pytest.raises(ValueError, match="three dimensions"):
            sr.transform_stack(bad, tmats=np.zeros((4, 3, 3)))


@pytest.mark.parametrize("shape", [(0, 20, 20), (20, 0, 20)])
def test_empty_stack(shape):
    sr = StackReg(StackReg.TRANSLATION)
    axis = shape.index(0)
    with pytest.raises(ValueError, match="no frames"):
        sr.register_stack(np.zeros(shape), axis=axis)
    with pytest.raises(ValueError, match="no frames"):
        sr.transform_stack(np.zeros(shape), axis=axis, tmats=np.zeros((0, 3, 3)))


def test_stack_with_empty_frames():
    with pytest.raises(ValueError, match="at least 3 x 3"):
        StackReg(StackReg.TRANSLATION).register_stack(
            np.zeros((4, 0, 20)), reference="first", suppress_axis_warning=True
        )


@pytest.mark.parametrize("axis", [3, -4, 10])
def test_axis_out_of_bounds(stack, axis):
    with pytest.raises(ValueError, match="out of bounds"):
        StackReg(StackReg.TRANSLATION).register_stack(stack, axis=axis)


def test_axis_must_be_integer(stack):
    with pytest.raises(TypeError, match="axis"):
        StackReg(StackReg.TRANSLATION).register_stack(stack, axis=0.0)


def test_unknown_reference(stack):
    with pytest.raises(ValueError, match="Unknown reference"):
        StackReg(StackReg.TRANSLATION).register_stack(stack, reference="last")


def test_bilinear_previous_is_refused(stack):
    with pytest.raises(ValueError, match="Bilinear"):
        StackReg(StackReg.BILINEAR).register_stack(stack, reference="previous")


@pytest.mark.parametrize(
    ("argument", "value"),
    [("n_frames", 0), ("n_frames", -1), ("moving_average", -1)],
)
def test_invalid_counts_are_rejected(stack, argument, value):
    with pytest.raises(ValueError, match=argument):
        StackReg(StackReg.TRANSLATION).register_stack(
            stack, reference="first", **{argument: value}
        )


def test_counts_must_be_integers(stack):
    with pytest.raises(TypeError, match="n_frames"):
        StackReg(StackReg.TRANSLATION).register_stack(
            stack, reference="first", n_frames=1.5
        )


def test_n_frames_larger_than_stack(stack):
    with pytest.raises(ValueError, match="exceeds"):
        StackReg(StackReg.TRANSLATION).register_stack(
            stack, reference="first", n_frames=5
        )


def test_tmats_must_match_stack(stack):
    sr = StackReg(StackReg.TRANSLATION)
    with pytest.raises(ValueError, match="does not match stack length"):
        sr.transform_stack(stack, tmats=np.repeat(np.identity(3)[None], 3, axis=0))
    with pytest.raises(ValueError, match=r"shape \(N, 3, 3\)"):
        sr.transform_stack(stack, tmats=np.zeros((4, 4, 4)))
    with pytest.raises(ValueError, match=r"shape \(N, 3, 3\)"):
        sr.transform_stack(stack, tmats=np.identity(3))


def test_nan_frame_in_stack_is_reported(stack):
    bad = stack.copy()
    bad[2, 3, 3] = np.nan
    with pytest.raises(ValueError, match="NaN or infinite"):
        StackReg(StackReg.TRANSLATION).register_stack(bad, reference="previous")


@pytest.mark.parametrize("n_jobs", [0, 1.5, "all", True])
def test_invalid_n_jobs(stack, n_jobs):
    error = ValueError if n_jobs == 0 else TypeError
    with pytest.raises(error, match="n_jobs"):
        StackReg(StackReg.TRANSLATION).register_stack(stack, n_jobs=n_jobs)


def test_invalid_out_argument(stack):
    sr = StackReg(StackReg.TRANSLATION)
    tmats = np.repeat(np.identity(3)[None], 4, axis=0)
    with pytest.raises(ValueError, match="shape"):
        sr.transform_stack(stack, tmats=tmats, out=np.empty((4, 32, 41)))
    with pytest.raises(TypeError, match="dtype"):
        sr.transform_stack(stack, tmats=tmats, out=np.empty(stack.shape, dtype=bool))
    with pytest.raises(TypeError, match="numpy array"):
        sr.transform_stack(stack, tmats=tmats, out=stack.tolist())
    readonly = np.empty(stack.shape)
    readonly.setflags(write=False)
    with pytest.raises(ValueError, match="read-only"):
        sr.transform_stack(stack, tmats=tmats, out=readonly)


def test_progress_callback_must_be_callable(stack):
    with pytest.raises(TypeError, match="callable"):
        StackReg(StackReg.TRANSLATION).register_stack(stack, progress_callback=3)
