"""Behaviour of the resampling step (``StackReg.transform``)."""

from __future__ import annotations

import numpy as np
import pytest

from pystackreg import StackReg
from pystackreg.util import to_int_dtype, to_uint16

from helpers import TRANSFORMATIONS, Scene


@pytest.fixture(scope="module")
def image():
    # single-precision values: TurboReg stores images as float
    return Scene((70, 90), seed=5).render().astype(np.float32).astype(np.float64)


def translation(dx, dy):
    return np.array([[1.0, 0.0, dx], [0.0, 1.0, dy], [0.0, 0.0, 1.0]])


@pytest.mark.parametrize("name", list(TRANSFORMATIONS))
def test_identity_reproduces_image(image, name):
    """Cubic B-spline interpolation is exact at the sample positions."""
    sr = StackReg(TRANSFORMATIONS[name])
    out = sr.transform(image, np.identity(4 if name == "bilinear" else 3))
    assert out.dtype == np.float64
    assert out.shape == image.shape
    np.testing.assert_allclose(out, image, rtol=2e-6)


def test_output_has_single_precision(image):
    """As in TurboReg, results are rounded to single precision."""
    out = StackReg(StackReg.AFFINE).transform(image, translation(0.3, 0.6))
    np.testing.assert_array_equal(out, out.astype(np.float32))


def test_integer_translation_moves_pixels(image):
    """A shift by whole pixels relocates samples: out[y, x] = img[y+dy, x+dx];
    positions without data are set to 0."""
    out = StackReg(StackReg.TRANSLATION).transform(image, translation(5.0, -3.0))
    np.testing.assert_allclose(out[3:, :-5], image[:-3, 5:], rtol=2e-6)
    np.testing.assert_array_equal(out[:3, :], 0.0)
    np.testing.assert_array_equal(out[:, -5:], 0.0)


def test_translation_ignores_linear_part(image):
    """A TRANSLATION object only uses the last column of the matrix."""
    tmat = np.array([[0.5, 0.3, 2.0], [0.1, 2.0, 1.0], [0.0, 0.0, 1.0]])
    np.testing.assert_array_equal(
        StackReg(StackReg.TRANSLATION).transform(image, tmat),
        StackReg(StackReg.TRANSLATION).transform(image, translation(2.0, 1.0)),
    )


def test_affine_types_share_the_resampling_code(image):
    tmat = np.array([[1.02, 0.03, -4.2], [-0.02, 0.97, 3.1], [0.0, 0.0, 1.0]])
    expected = StackReg(StackReg.AFFINE).transform(image, tmat)
    for code in (StackReg.RIGID_BODY, StackReg.SCALED_ROTATION):
        np.testing.assert_array_equal(StackReg(code).transform(image, tmat), expected)


def test_bilinear_with_zero_cross_term_equals_affine(image):
    affine = np.array([[1.02, 0.03, -4.2], [-0.02, 0.97, 3.1], [0.0, 0.0, 1.0]])
    bilinear = np.identity(4)
    bilinear[0] = [1.02, 0.03, 0.0, -4.2]
    bilinear[1] = [-0.02, 0.97, 0.0, 3.1]
    np.testing.assert_allclose(
        StackReg(StackReg.BILINEAR).transform(image, bilinear),
        StackReg(StackReg.AFFINE).transform(image, affine),
        rtol=1e-6,
        atol=1e-3,
    )


@pytest.mark.parametrize(
    "tmat",
    [
        translation(-2.25, 4.5),
        np.array([[1.02, 0.03, -4.2], [-0.02, 0.97, 3.1], [0.0, 0.0, 1.0]]),
        np.array([[0.9, -0.2, 9.0], [0.2, 0.9, -6.0], [0.0, 0.0, 1.0]]),
    ],
    ids=["translation", "affine", "rotation"],
)
def test_agrees_with_scipy_cubic_spline(image, tmat):
    """Independent implementation: SciPy's cubic spline interpolation with
    half-sample symmetric boundaries gives the same values wherever TurboReg
    produces data (nearest source pixel inside the image)."""
    ndimage = pytest.importorskip("scipy.ndimage")
    rows, cols = image.shape
    y, x = np.mgrid[0:rows, 0:cols].astype(float)
    src_x = tmat[0, 0] * x + tmat[0, 1] * y + tmat[0, 2]
    src_y = tmat[1, 0] * x + tmat[1, 1] * y + tmat[1, 2]
    expected = ndimage.map_coordinates(image, [src_y, src_x], order=3, mode="reflect")

    # TurboReg produces data where the nearest source pixel lies inside the
    # image, i.e. for source coordinates in [-0.5, n - 0.5). It accumulates
    # coordinates incrementally, so positions within rounding distance of
    # that boundary are not classified here.
    eps = 1e-6
    inside = (
        (src_x > -0.5 + eps)
        & (src_x < cols - 0.5 - eps)
        & (src_y > -0.5 + eps)
        & (src_y < rows - 0.5 - eps)
    )
    outside = (
        (src_x < -0.5 - eps)
        | (src_x > cols - 0.5 + eps)
        | (src_y < -0.5 - eps)
        | (src_y > rows - 0.5 + eps)
    )
    out = StackReg(StackReg.AFFINE).transform(image, tmat)
    assert inside.sum() > 1000
    assert outside.sum() > 100
    np.testing.assert_array_equal(out[outside], 0.0)
    np.testing.assert_allclose(out[inside], expected[inside], rtol=1e-6, atol=1e-3)


def test_interpolation_overshoots_at_edges():
    """Expected behaviour, not a defect (issues #8, #32): cubic interpolation
    of a step edge at sub-pixel positions under- and overshoots the range of
    the data, so a non-negative image can acquire negative values."""
    img = np.zeros((32, 32))
    img[:, 16:] = 1000.0
    out = StackReg(StackReg.TRANSLATION).transform(img, translation(0.5, 0.0))
    assert out.min() < -10
    assert out.max() > 1010
    # converting back to an integer type clips the excursions
    clipped = to_uint16(out)
    assert clipped.dtype == np.uint16
    assert clipped.min() == 0
    assert to_int_dtype(out, np.uint8).max() == 255


def test_mirror_boundary_at_image_border():
    """Samples up to half a pixel outside the image are extrapolated by
    mirroring; beyond that the output is 0."""
    img = np.full((16, 16), 500.0)
    out = StackReg(StackReg.TRANSLATION).transform(img, translation(-0.4, 0.0))
    np.testing.assert_allclose(out, 500.0, rtol=1e-6)
    out = StackReg(StackReg.TRANSLATION).transform(img, translation(-0.6, 0.0))
    np.testing.assert_array_equal(out[:, 0], 0.0)
    np.testing.assert_allclose(out[:, 1:], 500.0, rtol=1e-6)


@pytest.mark.parametrize("shape", [(1, 1), (1, 7), (7, 1), (2, 2), (3, 5)])
def test_small_images_can_be_transformed(shape):
    img = np.arange(np.prod(shape), dtype=float).reshape(shape) + 1.0
    out = StackReg(StackReg.AFFINE).transform(img, np.identity(3))
    np.testing.assert_allclose(out, img, rtol=1e-6)


def test_nan_spreads_over_the_whole_image(image):
    """The B-spline prefilter is recursive: one NaN contaminates every pixel."""
    img = image.copy()
    img[10, 10] = np.nan
    out = StackReg(StackReg.AFFINE).transform(img, translation(0.5, 0.5))
    assert np.isnan(out[out != 0]).all()


@pytest.mark.parametrize("dtype", ["uint8", "uint16", "int32", "float32", "float64"])
def test_input_dtype_does_not_matter(dtype):
    values = np.round(Scene((40, 50), seed=2).render() / 16.0)
    assert values.max() <= 255
    tmat = np.array([[1.01, 0.02, -1.2], [-0.02, 0.99, 0.7], [0.0, 0.0, 1.0]])
    sr = StackReg(StackReg.AFFINE)
    np.testing.assert_array_equal(
        sr.transform(values.astype(dtype), tmat), sr.transform(values, tmat)
    )


def test_input_layout_does_not_matter(image):
    tmat = np.array([[1.01, 0.02, -1.2], [-0.02, 0.99, 0.7], [0.0, 0.0, 1.0]])
    sr = StackReg(StackReg.AFFINE)
    expected = sr.transform(image, tmat)

    readonly = image.copy()
    readonly.setflags(write=False)
    strided = np.repeat(np.repeat(image, 2, axis=0), 3, axis=1)[::2, ::3]
    assert not strided.flags.c_contiguous
    for variant in (np.asfortranarray(image), readonly, strided, image.tolist()):
        np.testing.assert_array_equal(sr.transform(variant, tmat), expected)
    np.testing.assert_array_equal(readonly, image)  # input left untouched


def test_set_matrix_and_get_matrix_roundtrip(image):
    for name, code in TRANSFORMATIONS.items():
        size = 4 if name == "bilinear" else 3
        tmat = np.identity(size)
        tmat[0, -1], tmat[1, -1] = 1.5, -2.5
        sr = StackReg(code)
        assert not sr.is_registered()
        sr.set_matrix(tmat)
        assert sr.is_registered()
        np.testing.assert_array_equal(sr.get_matrix(), tmat)
        np.testing.assert_array_equal(sr.transform(image), sr.transform(image, tmat))
