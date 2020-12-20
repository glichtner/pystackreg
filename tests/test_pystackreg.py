import os
import pytest
import numpy as np
from pystackreg import StackReg
import tifffile as tif

BASE_PATH = "tests/data/"
prefix = "pc12"


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


@pytest.fixture
def stack_unregistered():
    return load_file("unregistered")


@pytest.fixture
def stack_translation():
    return load_file("translation")


@pytest.fixture
def stack_rigid_body():
    return load_file("rigid-body")


@pytest.fixture
def stack_scaled_roation():
    return load_file("scaled-rotation")


@pytest.fixture
def stack_affine():
    return load_file("affine")


@pytest.fixture
def stack_bilinear():
    return load_file("bilinear")


def load_file(transformation):

    if transformation == "unregistered":
        fname = prefix + "-unreg"
    else:
        fname = prefix + "-reg-" + transformation

    return tif.imread(os.path.join(BASE_PATH, fname + ".tif"))


@pytest.fixture(
    params=[
        StackReg.TRANSLATION,
        StackReg.RIGID_BODY,
        StackReg.SCALED_ROTATION,
        StackReg.AFFINE,
        StackReg.BILINEAR,
    ]
)
def stack(request):
    transformation_map = {
        StackReg.TRANSLATION: "translation",
        StackReg.RIGID_BODY: "rigid-body",
        StackReg.SCALED_ROTATION: "scaled-rotation",
        StackReg.AFFINE: "affine",
        StackReg.BILINEAR: "bilinear",
    }

    return {
        "transformation": request.param,
        "registered": load_file(transformation_map[request.param]),
    }


def test_registration_transformation(stack, stack_unregistered):
    sr = StackReg(stack["transformation"])
    reference = (
        "previous" if not stack["transformation"] == StackReg.BILINEAR else "first"
    )
    out = sr.register_transform_stack(stack_unregistered, reference=reference)

    np.testing.assert_allclose(
        to_uint16(out), to_uint16(stack["registered"]), rtol=1e-7, atol=1
    )


@pytest.mark.parametrize("frame_axis", [1, 2])
def test_different_axis(stack, stack_unregistered, frame_axis):
    stack["registered"] = np.moveaxis(stack["registered"], 0, frame_axis)
    stack_unregistered = np.moveaxis(stack_unregistered, 0, frame_axis)

    sr = StackReg(stack["transformation"])
    reference = (
        "previous" if not stack["transformation"] == StackReg.BILINEAR else "first"
    )
    out = sr.register_transform_stack(
        stack_unregistered, reference=reference, axis=frame_axis
    )

    assert out.shape == stack["registered"].shape
    np.testing.assert_allclose(
        to_uint16(out), to_uint16(stack["registered"]), rtol=1e-7, atol=1,
    )
