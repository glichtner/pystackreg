"""Coordinate transformation (issue #42) and transformation amplitude
(issue #41)."""

from __future__ import annotations

import numpy as np
import pytest

from pystackreg import StackReg
from pystackreg.util import max_displacement

from helpers import (
    TRANSFORMATIONS,
    Scene,
    known_transformations,
    rotation_about_center,
)

SHAPE = (96, 128)


def blob(shape, x0, y0, sigma=3.0):
    y, x = np.mgrid[0 : shape[0], 0 : shape[1]].astype(float)
    return 1000.0 * np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / (2 * sigma**2))


def centroid(img):
    img = np.clip(img, 0, None)
    y, x = np.mgrid[0 : img.shape[0], 0 : img.shape[1]]
    return np.array([(img * x).sum(), (img * y).sum()]) / img.sum()


@pytest.mark.parametrize("name", list(TRANSFORMATIONS))
def test_points_move_like_image_content(name):
    """A feature at ``p`` in the moving image is found at
    ``transform_points(p)`` in the transformed image."""
    tmat = known_transformations(SHAPE)[name]
    position = np.array([70.0, 40.0])  # (x, y)
    mov = blob(SHAPE, *position)
    sr = StackReg(TRANSFORMATIONS[name])
    out = sr.transform(mov, tmat)
    predicted = sr.transform_points(position, tmat)
    np.testing.assert_allclose(centroid(out), predicted, atol=0.05)
    # the plain matrix product is the opposite direction
    back = sr.transform_points(predicted, tmat, inverse=True)
    np.testing.assert_allclose(back, position, atol=1e-8)


@pytest.mark.parametrize("name", ["rigid_body", "scaled_rotation", "affine"])
def test_inverse_direction_is_matrix_product(name):
    tmat = known_transformations(SHAPE)[name]
    pts = np.array([[0.0, 0.0], [10.0, 20.0], [127.0, 95.0]])
    result = StackReg(TRANSFORMATIONS[name]).transform_points(pts, tmat, inverse=True)
    homogeneous = np.column_stack([pts, np.ones(len(pts))])
    np.testing.assert_allclose(result, (tmat @ homogeneous.T).T[:, :2], atol=1e-12)


@pytest.mark.parametrize("name", list(TRANSFORMATIONS))
def test_roundtrip_and_shapes(name):
    tmat = known_transformations(SHAPE)[name]
    sr = StackReg(TRANSFORMATIONS[name])
    grid = np.stack(np.meshgrid(np.linspace(0, 127, 7), np.linspace(0, 95, 5)), axis=-1)
    assert grid.shape == (5, 7, 2)
    forward = sr.transform_points(grid, tmat)
    assert forward.shape == grid.shape
    np.testing.assert_allclose(
        sr.transform_points(forward, tmat, inverse=True), grid, atol=1e-7
    )
    single = sr.transform_points([3.0, 4.0], tmat)
    assert single.shape == (2,)
    assert sr.transform_points(np.empty((0, 2)), tmat).shape == (0, 2)


def test_uses_last_registration_by_default():
    scene = Scene(SHAPE)
    truth = known_transformations(SHAPE)["affine"]
    sr = StackReg(StackReg.AFFINE)
    tmat = sr.register(scene.render(), scene.render(truth))
    pts = np.array([[20.0, 30.0], [100.0, 60.0]])
    np.testing.assert_array_equal(
        sr.transform_points(pts), sr.transform_points(pts, tmat)
    )
    # ground truth: reference position -> position in the moving image
    homogeneous = np.column_stack([pts, np.ones(2)])
    expected = (truth @ homogeneous.T).T[:, :2]
    np.testing.assert_allclose(
        sr.transform_points(pts, inverse=True), expected, atol=0.01
    )


def test_translation_uses_only_the_shift():
    tmat = np.array([[0.5, 0.2, 3.0], [0.1, 2.0, -4.0], [0.0, 0.0, 1.0]])
    result = StackReg(StackReg.TRANSLATION).transform_points([[10.0, 10.0]], tmat)
    np.testing.assert_allclose(result, [[7.0, 14.0]])


def test_row_column_order_conversion():
    """Points from NumPy indexing are (row, column) and must be flipped."""
    tmat = np.array([[1.0, 0.0, 5.0], [0.0, 1.0, -3.0], [0.0, 0.0, 1.0]])
    sr = StackReg(StackReg.TRANSLATION)
    img = np.zeros((40, 50))
    img[20, 30] = 1.0
    out = sr.transform(img, tmat)
    row_col = np.argwhere(img == 1.0).astype(float)
    moved = sr.transform_points(row_col[:, ::-1], tmat)[:, ::-1]
    np.testing.assert_array_equal(np.argwhere(out > 0.5), moved)


def test_invalid_points():
    sr = StackReg(StackReg.AFFINE)
    for bad in (np.zeros((4, 3)), 1.0, np.zeros(3)):
        with pytest.raises(ValueError, match=r"\(\.\.\., 2\)"):
            sr.transform_points(bad, np.identity(3))


def test_singular_bilinear_transformation_is_reported():
    tmat = np.identity(4)
    tmat[0] = [0.0, 0.0, 0.0, 5.0]  # maps everything onto x = 5
    with pytest.raises(RuntimeError, match="invert the bilinear"):
        StackReg(StackReg.BILINEAR).transform_points([[1.0, 2.0]], tmat)


# -- max_displacement ------------------------------------------------------------


def test_max_displacement_identity_and_translation():
    assert max_displacement(np.identity(3), SHAPE) == 0.0
    assert max_displacement(np.identity(4), SHAPE) == 0.0
    tmat = np.array([[1.0, 0.0, 3.0], [0.0, 1.0, -4.0], [0.0, 0.0, 1.0]])
    assert max_displacement(tmat, SHAPE) == pytest.approx(5.0)


def test_max_displacement_rotation():
    """Rotation about the centre: the corners move by 2 r sin(angle / 2)."""
    angle = 10.0
    tmat = rotation_about_center(angle, SHAPE)
    radius = 0.5 * np.hypot(SHAPE[0] - 1, SHAPE[1] - 1)
    expected = 2 * radius * np.sin(np.deg2rad(angle) / 2)
    assert max_displacement(tmat, SHAPE) == pytest.approx(expected)


def test_max_displacement_is_attained_at_a_corner():
    rng = np.random.default_rng(1)
    y, x = np.mgrid[0 : SHAPE[0], 0 : SHAPE[1]].astype(float)
    for name in ("affine", "bilinear"):
        tmat = known_transformations(SHAPE)[name].copy()
        tmat[:2] += rng.normal(0, 1e-3, tmat[:2].shape)
        if name == "affine":
            u = tmat[0, 0] * x + tmat[0, 1] * y + tmat[0, 2]
            v = tmat[1, 0] * x + tmat[1, 1] * y + tmat[1, 2]
        else:
            u = tmat[0, 0] * x + tmat[0, 1] * y + tmat[0, 2] * x * y + tmat[0, 3]
            v = tmat[1, 0] * x + tmat[1, 1] * y + tmat[1, 2] * x * y + tmat[1, 3]
        brute_force = np.hypot(u - x, v - y).max()
        assert max_displacement(tmat, SHAPE) == pytest.approx(brute_force)


def test_max_displacement_for_a_stack_of_matrices():
    tmats = np.repeat(np.identity(3)[None], 4, axis=0)
    tmats[:, 0, 2] = [0.0, 1.0, 2.0, 30.0]
    result = max_displacement(tmats, SHAPE)
    np.testing.assert_allclose(result, [0.0, 1.0, 2.0, 30.0])
    assert list(np.flatnonzero(result > 25)) == [3]  # the implausible frame


def test_max_displacement_rejects_other_shapes():
    with pytest.raises(ValueError, match="shape"):
        max_displacement(np.identity(2), SHAPE)
    with pytest.raises(ValueError, match="shape"):
        max_displacement(np.zeros((2, 2, 3, 3)), SHAPE)
