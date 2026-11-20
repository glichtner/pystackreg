"""Numerical regression tests against pystackreg 0.2.8.

``tests/data/baseline_v0.2.8.npz`` holds results of the unmodified 0.2.8
release (see ``generate_baseline.py``). The cases are recomputed with the
current build and compared.

On the platform the baseline was generated on (Linux x86-64, GCC) the current
code reproduces it bit for bit; set ``PYSTACKREG_EXACT_BASELINE=1`` to enforce
that. Elsewhere, differences in the math library (``atan2``, ``sin``, ``cos``,
``exp``) can change the last bits, and the tests use tight tolerances instead.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from pystackreg import StackReg

from generate_baseline import (
    DTYPES,
    EXPLICIT_MATRICES,
    PAIR_SHAPES,
    as_dtype,
)
from helpers import TRANSFORMATIONS, Scene, drifting_stack, known_transformations

EXACT = os.environ.get("PYSTACKREG_EXACT_BASELINE") == "1"

#: absolute tolerance for matrix entries and landmark coordinates (pixels)
MATRIX_ATOL = 1e-6
#: relative tolerance for intensities of transformed images
IMAGE_RTOL = 1e-5


def assert_matches(actual, expected, *, image=False):
    if EXACT:
        np.testing.assert_array_equal(actual, expected)
    elif image:
        scale = float(np.abs(expected).max())
        np.testing.assert_allclose(actual, expected, rtol=0, atol=IMAGE_RTOL * scale)
    else:
        np.testing.assert_allclose(actual, expected, rtol=0, atol=MATRIX_ATOL)


def test_baseline_was_generated_with_0_2_8(baseline):
    assert str(baseline["meta/pystackreg_version"]) == "0.2.8"


@pytest.mark.parametrize("shape", PAIR_SHAPES, ids=lambda s: f"{s[0]}x{s[1]}")
@pytest.mark.parametrize("name", list(TRANSFORMATIONS))
def test_pair_registration(baseline, shape, name):
    scene = Scene(shape)
    ref = scene.render()
    mov = scene.render(known_transformations(shape)[name])
    key = f"pair/{shape[0]}x{shape[1]}/{name}"

    sr = StackReg(TRANSFORMATIONS[name])
    assert_matches(sr.register(ref, mov), baseline[key + "/matrix"])
    refpts, movpts = sr.get_points()
    assert_matches(refpts, baseline[key + "/refpts"])
    assert_matches(movpts, baseline[key + "/movpts"])
    if key + "/transformed" in baseline:
        assert_matches(sr.transform(mov), baseline[key + "/transformed"], image=True)


@pytest.mark.parametrize("case", list(EXPLICIT_MATRICES))
def test_transform_with_explicit_matrix(baseline, case):
    name, tmat = EXPLICIT_MATRICES[case]
    img = Scene((70, 90), seed=5).render()
    out = StackReg(TRANSFORMATIONS[name]).transform(img, tmat)
    assert_matches(out, baseline["explicit/" + case], image=True)


@pytest.mark.parametrize("reference", ["previous", "first", "mean"])
@pytest.mark.parametrize("name", list(TRANSFORMATIONS))
def test_stack_registration(baseline, pc12, name, reference):
    if name == "bilinear" and reference == "previous":
        pytest.skip("bilinear transformations cannot be chained")
    sr = StackReg(TRANSFORMATIONS[name])
    tmats = sr.register_stack(pc12, reference=reference)
    assert_matches(tmats, baseline[f"pc12/{name}/{reference}"])


def test_stack_registration_first_n_frames(baseline, pc12):
    tmats = StackReg(StackReg.RIGID_BODY).register_stack(
        pc12, reference="first", n_frames=3
    )
    assert_matches(tmats, baseline["pc12/rigid_body/first_n3"])


@pytest.mark.parametrize("reference", ["first", "mean"])
def test_stack_registration_moving_average(baseline, pc12, reference):
    # integer stack: the frame-wise moving average is exactly the one that
    # 0.2.8 computed from cumulative sums
    tmats = StackReg(StackReg.RIGID_BODY).register_stack(
        pc12, reference=reference, moving_average=3
    )
    assert_matches(tmats, baseline[f"pc12/rigid_body/{reference}_ma3"])


@pytest.mark.parametrize("dtype", DTYPES)
def test_stack_registration_dtypes(baseline, dtype):
    stack, _ = drifting_stack()
    tmats = StackReg(StackReg.RIGID_BODY).register_stack(
        as_dtype(stack, dtype), reference="previous"
    )
    assert_matches(tmats, baseline["dtype/" + dtype])
