"""Scientific validation: does the registration recover known transformations?

All tests use analytic scenes (``helpers.Scene``) for which a transformed copy
is computed exactly, so the expected transformation matrix is known without
any resampling error. Accuracy is expressed as the largest distance, in
pixels, between the positions to which the estimated and the true
transformation map the four image corners.
"""

from __future__ import annotations

import numpy as np
import pytest

from pystackreg import StackReg

from helpers import (
    TRANSFORMATIONS,
    Scene,
    corner_error,
    drifting_stack,
    known_transformations,
    rotation_about_center,
)

SHAPES = [(96, 128), (128, 96), (150, 110), (64, 64)]
AFFINE_KINDS = ["translation", "rigid_body", "scaled_rotation", "affine"]


@pytest.mark.parametrize("shape", SHAPES, ids=lambda s: f"{s[0]}x{s[1]}")
@pytest.mark.parametrize("name", list(TRANSFORMATIONS))
def test_recovers_known_transformation(name, shape):
    """Each transformation type is recovered to a few thousandths of a pixel
    (measured: < 0.001 px; the non-linear bilinear model is slightly less
    exact on the smallest image)."""
    scene = Scene(shape)
    truth = known_transformations(shape)[name]
    tmat = StackReg(TRANSFORMATIONS[name]).register(scene.render(), scene.render(truth))
    assert tmat.shape == truth.shape
    assert corner_error(tmat, truth, shape) < 0.01


def _smooth_scene(shape, dx=0.0, dy=0.0):
    """Smooth structure that is well sampled on a grid of a few pixels, with
    the structure displaced by ``(dx, dy)``."""
    height, width = shape
    y, x = np.mgrid[0:height, 0:width].astype(float)
    x, y = x - dx, y - dy
    size = min(height, width)
    blob = np.exp(
        -((x - width / 2) ** 2 + (y - height / 2) ** 2) / (2 * (size / 5) ** 2)
    )
    return 100 + 60 * blob + 20 * np.sin(0.5 * x) * np.cos(0.4 * y)


@pytest.mark.parametrize(
    "shape",
    [(16, 16), (12, 20), (20, 12), (23, 23), (23, 40)],
    ids=lambda s: f"{s[0]}x{s[1]}",
)
@pytest.mark.parametrize("name", AFFINE_KINDS)
def test_images_without_pyramid_level_are_registered(name, shape):
    """Images with a side shorter than 24 px have no pyramid level and are
    registered at full resolution only. Up to 0.2.8 the scaled-rotation and
    affine models read a value that was never set in this case; the result
    depended on the contents of memory and was usually the identity, which is
    also what TurboReg returns. The tolerance reflects how little information
    such images hold (measured: < 0.2 px; the identity is 0.94 px off)."""
    truth = np.identity(3)
    truth[0, 2], truth[1, 2] = 0.8, -0.5
    tmat = StackReg(TRANSFORMATIONS[name]).register(
        _smooth_scene(shape), _smooth_scene(shape, 0.8, -0.5)
    )
    assert corner_error(tmat, truth, shape) < 0.3


@pytest.mark.parametrize("name", list(TRANSFORMATIONS))
def test_identical_images_give_identity(name):
    ref = Scene((96, 128)).render()
    tmat = StackReg(TRANSFORMATIONS[name]).register(ref, ref)
    np.testing.assert_allclose(tmat, np.identity(tmat.shape[0]), atol=1e-6)


@pytest.mark.parametrize("shift", [0.05, 0.1, 0.25, 0.5, 0.75])
@pytest.mark.parametrize("name", ["translation", "rigid_body", "affine"])
def test_subpixel_translation(name, shift):
    scene = Scene((96, 128))
    truth = np.identity(3)
    truth[0, 2], truth[1, 2] = shift, -shift / 2
    tmat = StackReg(TRANSFORMATIONS[name]).register(scene.render(), scene.render(truth))
    assert corner_error(tmat, truth, scene.shape) < 1e-3


def test_matrix_maps_reference_to_moving_coordinates():
    """Sign and axis conventions of the returned matrix.

    The moving image shows the scene displaced by +3 px along x (columns) and
    -2 px along y (rows). The matrix acts on (x, y, 1) and maps a reference
    position onto the position of the same structure in the moving image, so
    its translation column is (+3, -2).
    """
    scene = Scene((96, 128))
    y, x = np.mgrid[0:96, 0:128].astype(float)
    ref = scene.evaluate(x, y)
    mov = scene.evaluate(x - 3.0, y + 2.0)  # structure moved to (x + 3, y - 2)
    tmat = StackReg(StackReg.TRANSLATION).register(ref, mov)
    np.testing.assert_allclose(tmat[:2, 2], [3.0, -2.0], atol=1e-3)
    np.testing.assert_array_equal(tmat[:2, :2], np.identity(2))
    np.testing.assert_array_equal(tmat[2], [0.0, 0.0, 1.0])


@pytest.mark.parametrize("name", list(TRANSFORMATIONS))
def test_transformed_image_matches_reference(name):
    """After register + transform the moving image coincides with the
    reference (away from the border, where data are missing)."""
    scene = Scene((96, 128))
    ref = scene.render()
    mov = scene.render(known_transformations(scene.shape)[name])
    out = StackReg(TRANSFORMATIONS[name]).register_transform(ref, mov)
    inner = (slice(12, -12), slice(12, -12))
    before = np.abs(mov - ref)[inner].max()
    after = np.abs(out - ref)[inner].max()
    assert before > 50
    # residual: interpolation error of the cubic spline on this scene
    assert after < 2.5
    assert after < before / 50


@pytest.mark.parametrize("name", AFFINE_KINDS)
def test_inverse_consistency(name):
    """Swapping reference and moving image yields the inverse transformation."""
    scene = Scene((96, 128))
    ref = scene.render()
    mov = scene.render(known_transformations(scene.shape)[name])
    sr = StackReg(TRANSFORMATIONS[name])
    forward = sr.register(ref, mov)
    backward = sr.register(mov, ref)
    assert corner_error(forward @ backward, np.identity(3), scene.shape) < 0.01


def test_matrix_structure_per_transformation_type():
    scene = Scene((96, 128))
    ref = scene.render()
    mov = scene.render(known_transformations(scene.shape)["affine"])

    tra = StackReg(StackReg.TRANSLATION).register(ref, mov)
    np.testing.assert_array_equal(tra[:2, :2], np.identity(2))

    rig = StackReg(StackReg.RIGID_BODY).register(ref, mov)[:2, :2]
    np.testing.assert_allclose(rig.T @ rig, np.identity(2), atol=1e-12)
    assert np.linalg.det(rig) == pytest.approx(1.0, abs=1e-12)

    sca = StackReg(StackReg.SCALED_ROTATION).register(ref, mov)[:2, :2]
    np.testing.assert_allclose(sca[0, 0], sca[1, 1], atol=1e-9)
    np.testing.assert_allclose(sca[0, 1], -sca[1, 0], atol=1e-9)

    for code in (StackReg.RIGID_BODY, StackReg.SCALED_ROTATION, StackReg.AFFINE):
        np.testing.assert_array_equal(
            StackReg(code).register(ref, mov)[2], [0.0, 0.0, 1.0]
        )
    bil = StackReg(StackReg.BILINEAR).register(ref, mov)
    assert bil.shape == (4, 4)
    np.testing.assert_array_equal(bil[2:], np.identity(4)[2:])


@pytest.mark.parametrize("name", AFFINE_KINDS)
def test_invariance_to_intensity_scaling_and_offset(name):
    """The estimate does not depend on the intensity scale or offset."""
    scene = Scene((96, 128))
    ref = scene.render()
    mov = scene.render(known_transformations(scene.shape)[name])
    sr = StackReg(TRANSFORMATIONS[name])
    expected = sr.register(ref, mov)
    np.testing.assert_allclose(sr.register(ref * 1e-3, mov * 1e-3), expected, atol=1e-5)
    np.testing.assert_allclose(
        sr.register(ref * 1e3 + 5e4, mov * 1e3 + 5e4), expected, atol=1e-5
    )


def test_last_row_and_column_are_ignored():
    """As in StackReg, registration does not use the last row and column."""
    scene = Scene((96, 128))
    ref = scene.render()
    mov = scene.render(known_transformations(scene.shape)["rigid_body"])
    sr = StackReg(StackReg.RIGID_BODY)
    expected = sr.register(ref, mov)
    ref2, mov2 = ref.copy(), mov.copy()
    ref2[-1, :] = 1e6
    mov2[:, -1] = -1e6
    np.testing.assert_array_equal(sr.register(ref2, mov2), expected)


@pytest.mark.parametrize("reference", ["previous", "first", "mean"])
@pytest.mark.parametrize("kind", ["translation", "rigid_body"])
def test_stack_drift_is_recovered(kind, reference):
    """Per-frame matrices of a drifting stack agree with the ground truth.

    For ``previous`` this validates the chaining of pairwise transformations:
    ``tmats[i]`` maps frame-0 coordinates onto frame-i coordinates.
    """
    stack, truth = drifting_stack(n_frames=8, kind=kind)
    tmats = StackReg(TRANSFORMATIONS[kind]).register_stack(stack, reference=reference)
    if reference == "mean":
        # the mean image defines its own coordinate frame: compare the
        # transformations relative to frame 0
        tmats = np.array([t @ np.linalg.inv(tmats[0]) for t in tmats])
    errors = [
        corner_error(t, e, stack.shape[1:]) for t, e in zip(tmats, truth, strict=True)
    ]
    assert max(errors) < 0.05


@pytest.mark.parametrize("reference", ["previous", "first", "mean"])
def test_registered_stack_is_aligned(reference):
    stack, _ = drifting_stack(n_frames=8, kind="rigid_body")
    out = StackReg(StackReg.RIGID_BODY).register_transform_stack(
        stack, reference=reference
    )
    inner = (slice(None), slice(10, -10), slice(10, -10))
    spread_before = np.ptp(stack[inner], axis=0).max()
    spread_after = np.ptp(out[inner], axis=0).max()
    assert spread_after < spread_before / 30


@pytest.mark.slow
@pytest.mark.parametrize("name", ["rigid_body", "affine"])
def test_large_image(name):
    """A large, non-square image (six pyramid levels)."""
    shape = (768, 1152)
    scene = Scene(shape, n_blobs=120, seed=11)
    truth = known_transformations(shape)[name]
    mov = scene.render(truth)
    sr = StackReg(TRANSFORMATIONS[name])
    tmat = sr.register(scene.render(), mov)
    assert corner_error(tmat, truth, shape) < 0.01
    out = sr.transform(mov)
    assert out.shape == shape
    assert np.isfinite(out).all()


def test_capture_range_is_limited():
    """Documented limitation: like TurboReg, the optimiser is local.

    Moderate displacements are recovered; a rotation by 90 degrees is not
    (see issue #30). The test pins both sides of that statement.
    """
    shape = (128, 128)
    scene = Scene(shape, seed=4)
    ref = scene.render()
    sr = StackReg(StackReg.RIGID_BODY)

    moderate = rotation_about_center(20.0, shape)
    assert (
        corner_error(sr.register(ref, scene.render(moderate)), moderate, shape) < 0.01
    )

    large = rotation_about_center(90.0, shape)
    assert corner_error(sr.register(ref, scene.render(large)), large, shape) > 10
