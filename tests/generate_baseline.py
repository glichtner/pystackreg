"""Generate the numerical regression baseline used by ``test_regression.py``.

The committed baseline (``tests/data/baseline_v0.2.8.npz``) was produced with
the unmodified pystackreg 0.2.8 release (git tag ``v0.2.8``), so that the tests
detect any unintended numerical change introduced afterwards. It was generated
on Linux x86-64 (GCC, CPython 3.12, NumPy 2.5).

Only re-run this script against a *trusted* build, and only when a deliberate,
documented numerical change makes a new baseline necessary::

    python tests/generate_baseline.py tests/data/baseline_vX.Y.Z.npz
"""

from __future__ import annotations

import os
import sys

import numpy as np
import tifffile

import pystackreg
from pystackreg import StackReg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from helpers import (
    TRANSFORMATIONS,
    Scene,
    drifting_stack,
    known_transformations,
)

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

#: image sizes (height, width) covered by the single-pair cases
PAIR_SHAPES = [(96, 128), (61, 47), (25, 33)]

#: explicit matrices applied with ``transform`` (no registration involved)
EXPLICIT_MATRICES = {
    "translation_integer": (
        "translation",
        np.array([[1.0, 0.0, 5.0], [0.0, 1.0, -3.0], [0.0, 0.0, 1.0]]),
    ),
    "translation_subpixel": (
        "translation",
        np.array([[1.0, 0.0, -2.25], [0.0, 1.0, 4.5], [0.0, 0.0, 1.0]]),
    ),
    "affine": (
        "affine",
        np.array([[1.02, 0.03, -4.2], [-0.02, 0.97, 3.1], [0.0, 0.0, 1.0]]),
    ),
    "bilinear": (
        "bilinear",
        np.array(
            [
                [1.01, 0.02, 1.0e-4, -1.5],
                [-0.01, 0.98, -5.0e-5, 2.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        ),
    ),
}

DTYPES = ["uint8", "uint16", "int32", "float32", "float64"]


def as_dtype(stack: np.ndarray, dtype: str) -> np.ndarray:
    """Convert the float test stack to ``dtype`` (scaled into uint8 range)."""
    if dtype == "uint8":
        return np.round(stack / stack.max() * 255).astype(np.uint8)
    if np.issubdtype(np.dtype(dtype), np.integer):
        return np.round(stack).astype(dtype)
    return stack.astype(dtype)


def generate() -> dict:
    out = {}

    # -- single image pairs with known transformations ----------------------
    for shape in PAIR_SHAPES:
        scene = Scene(shape)
        ref = scene.render()
        for name, tmat in known_transformations(shape).items():
            key = f"pair/{shape[0]}x{shape[1]}/{name}"
            mov = scene.render(tmat)
            sr = StackReg(TRANSFORMATIONS[name])
            out[key + "/matrix"] = sr.register(ref, mov)
            refpts, movpts = sr.get_points()
            out[key + "/refpts"] = refpts
            out[key + "/movpts"] = movpts
            if shape == PAIR_SHAPES[0]:
                # transform output is single precision by construction
                out[key + "/transformed"] = sr.transform(mov).astype(np.float32)

    # -- explicit matrices ---------------------------------------------------
    img = Scene((70, 90), seed=5).render()
    for key, (name, tmat) in EXPLICIT_MATRICES.items():
        sr = StackReg(TRANSFORMATIONS[name])
        out["explicit/" + key] = sr.transform(img, tmat).astype(np.float32)

    # -- real data: every transformation and reference mode -----------------
    pc12 = tifffile.imread(os.path.join(DATA, "pc12-unreg.tif"))
    for name, code in TRANSFORMATIONS.items():
        for reference in ("previous", "first", "mean"):
            if name == "bilinear" and reference == "previous":
                continue
            sr = StackReg(code)
            out[f"pc12/{name}/{reference}"] = sr.register_stack(
                pc12, reference=reference
            )
    sr = StackReg(StackReg.RIGID_BODY)
    out["pc12/rigid_body/first_n3"] = sr.register_stack(
        pc12, reference="first", n_frames=3
    )
    for reference in ("first", "mean"):
        sr = StackReg(StackReg.RIGID_BODY)
        out[f"pc12/rigid_body/{reference}_ma3"] = sr.register_stack(
            pc12, reference=reference, moving_average=3
        )
    # -- dtypes ---------------------------------------------------------------
    stack, _ = drifting_stack()
    for dtype in DTYPES:
        sr = StackReg(StackReg.RIGID_BODY)
        out["dtype/" + dtype] = sr.register_stack(
            as_dtype(stack, dtype), reference="previous"
        )

    return out


if __name__ == "__main__":
    target = sys.argv[1]
    result = generate()
    result["meta/pystackreg_version"] = np.array(pystackreg.__version__)
    result["meta/numpy_version"] = np.array(np.__version__)
    np.savez_compressed(target, **result)
    print(
        f"wrote {len(result)} arrays to {target} using pystackreg "
        f"{pystackreg.__version__} from {pystackreg.__file__}"
    )
