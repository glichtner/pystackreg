"""Comparison with images registered by the original ImageJ plugins.

The reference stacks in ``tests/data`` (``pc12-reg-*.tif``) were produced with
StackReg/TurboReg in ImageJ from ``pc12-unreg.tif`` using the macros in
``examples/data``. ImageJ stores the result as 16-bit integers, so agreement
is checked after rounding, with a tolerance of one grey level.
"""

from __future__ import annotations

import numpy as np
import pytest

from pystackreg import StackReg

from conftest import DATA_DIR

tifffile = pytest.importorskip("tifffile")

REFERENCE_FILES = {
    StackReg.TRANSLATION: "translation",
    StackReg.RIGID_BODY: "rigid-body",
    StackReg.SCALED_ROTATION: "scaled-rotation",
    StackReg.AFFINE: "affine",
    StackReg.BILINEAR: "bilinear",
}


def np_school_round(x):
    """
    Rounding function that always round .5 value up (as opposed to standard
    python/numpy behaviour, which round .5 values to the next even number.
    This function here mimics the behavior of ImageJ for comparability of
    StackReg/TurboReg registrations from ImageJ and pystackreg.

    :param x: (float) number
    :return: (float) rounded number (.5 always rounded up)
    """
    x = np.array(x)
    r = x.copy()
    idx = x - np.floor(x) < 0.5
    r[idx] = np.floor(x[idx])
    r[~idx] = np.ceil(x[~idx])

    return r


def to_uint16(img):
    return np_school_round(img.clip(min=0, max=65535)).astype(np.uint16)


@pytest.fixture(params=list(REFERENCE_FILES), ids=list(REFERENCE_FILES.values()))
def stack(request):
    name = REFERENCE_FILES[request.param]
    return {
        "transformation": request.param,
        "registered": tifffile.imread(DATA_DIR / f"pc12-reg-{name}.tif"),
    }


def reference_mode(transformation):
    # StackReg aligns every slice to the previous one. It does not offer the
    # bilinear transformation; that reference was made with TurboReg against
    # the first slice.
    return "first" if transformation == StackReg.BILINEAR else "previous"


def test_registration_transformation(stack, pc12):
    sr = StackReg(stack["transformation"])
    out = sr.register_transform_stack(
        pc12, reference=reference_mode(stack["transformation"])
    )

    np.testing.assert_allclose(
        to_uint16(out), to_uint16(stack["registered"]), rtol=1e-7, atol=1
    )


@pytest.mark.parametrize("frame_axis", [1, 2])
def test_different_axis(stack, pc12, frame_axis):
    registered = np.moveaxis(stack["registered"], 0, frame_axis)
    unregistered = np.moveaxis(pc12, 0, frame_axis)

    sr = StackReg(stack["transformation"])
    out = sr.register_transform_stack(
        unregistered,
        reference=reference_mode(stack["transformation"]),
        axis=frame_axis,
    )

    assert out.shape == registered.shape
    np.testing.assert_allclose(
        to_uint16(out),
        to_uint16(registered),
        rtol=1e-7,
        atol=1,
    )


@pytest.mark.parametrize("n_jobs", [2, -1])
def test_parallel_registration_matches_imagej(stack, pc12, n_jobs):
    sr = StackReg(stack["transformation"])
    out = sr.register_transform_stack(
        pc12, reference=reference_mode(stack["transformation"]), n_jobs=n_jobs
    )
    np.testing.assert_allclose(
        to_uint16(out), to_uint16(stack["registered"]), rtol=1e-7, atol=1
    )
