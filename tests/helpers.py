"""Deterministic synthetic data shared by the test suite and the baseline
generator.

The scene is an analytic function (a sum of Gaussian blobs on a smooth
background), so a geometrically transformed copy can be produced by evaluating
the function at transformed coordinates instead of by resampling. This gives an
exact ground truth for the transformation that a registration should recover.

No random number generator is used: blob parameters come from a fixed linear
congruential sequence, so the data are identical on every platform and NumPy
version (up to the last-bit behaviour of ``exp``).
"""

from __future__ import annotations

import numpy as np

from pystackreg import StackReg

TRANSFORMATIONS = {
    "translation": StackReg.TRANSLATION,
    "rigid_body": StackReg.RIGID_BODY,
    "scaled_rotation": StackReg.SCALED_ROTATION,
    "affine": StackReg.AFFINE,
    "bilinear": StackReg.BILINEAR,
}


def _lcg(seed: int, n: int) -> np.ndarray:
    """Return ``n`` reproducible pseudo-random numbers in [0, 1)."""
    out = np.empty(n)
    state = (seed * 2654435761 + 12345) % 2**32
    for i in range(n):
        state = (1664525 * state + 1013904223) % 2**32
        out[i] = state / 2**32
    return out


class Scene:
    """Analytic test scene with support for exact geometric transformation.

    :param shape: ``(height, width)`` of the rendered images
    :param n_blobs: number of Gaussian blobs
    :param seed: seed of the deterministic parameter sequence
    """

    def __init__(self, shape=(96, 128), n_blobs=40, seed=1):
        self.shape = tuple(shape)
        height, width = self.shape
        u = _lcg(seed, 4 * n_blobs).reshape(4, n_blobs)
        # keep blobs away from the border so that moderate transformations do
        # not move a large part of the structure out of the field of view
        self._cx = (0.15 + 0.7 * u[0]) * (width - 1)
        self._cy = (0.15 + 0.7 * u[1]) * (height - 1)
        self._sigma = (0.02 + 0.04 * u[2]) * min(height, width)
        self._amp = 200.0 + 800.0 * u[3]

    def evaluate(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Evaluate the scene at coordinates ``x`` (column) and ``y`` (row)."""
        height, width = self.shape
        out = 100.0 + 20.0 * x / width + 10.0 * y / height
        for cx, cy, sigma, amp in zip(self._cx, self._cy, self._sigma, self._amp):
            out = out + amp * np.exp(
                -((x - cx) ** 2 + (y - cy) ** 2) / (2.0 * sigma**2)
            )
        return out

    def render(self, tmat: np.ndarray | None = None) -> np.ndarray:
        """Render the scene, optionally seen through a transformation.

        ``tmat`` uses the pystackreg convention: it maps coordinates of the
        reference image onto coordinates of the moving image. The returned
        image is therefore the *moving* image for which a registration against
        ``render()`` should return ``tmat``.

        :param tmat: 3x3 matrix acting on ``(x, y, 1)`` or 4x4 bilinear matrix
            acting on ``(x, y, x*y, 1)``; ``None`` renders the reference
        """
        height, width = self.shape
        y, x = np.mgrid[0:height, 0:width].astype(float)
        if tmat is None:
            return self.evaluate(x, y)
        inv = inverse_map(np.asarray(tmat, dtype=float), x, y)
        return self.evaluate(*inv)


def inverse_map(tmat: np.ndarray, x: np.ndarray, y: np.ndarray):
    """Map moving-image coordinates back to reference coordinates.

    For 3x3 matrices this is the matrix inverse. The bilinear map has no
    closed-form matrix inverse and is inverted by Newton iteration.
    """
    if tmat.shape == (3, 3):
        inv = np.linalg.inv(tmat)
        return (
            inv[0, 0] * x + inv[0, 1] * y + inv[0, 2],
            inv[1, 0] * x + inv[1, 1] * y + inv[1, 2],
        )
    a, b = tmat[0], tmat[1]
    u, v = x.copy(), y.copy()
    for _ in range(50):
        fu = a[0] * u + a[1] * v + a[2] * u * v + a[3] - x
        fv = b[0] * u + b[1] * v + b[2] * u * v + b[3] - y
        j11, j12 = a[0] + a[2] * v, a[1] + a[2] * u
        j21, j22 = b[0] + b[2] * v, b[1] + b[2] * u
        det = j11 * j22 - j12 * j21
        du = (j22 * fu - j12 * fv) / det
        dv = (j11 * fv - j21 * fu) / det
        u, v = u - du, v - dv
        if max(np.abs(du).max(), np.abs(dv).max()) < 1e-12:
            break
    return u, v


def rotation_about_center(angle_deg, shape, scale=1.0, shift=(0.0, 0.0)):
    """3x3 matrix: rotation and isotropic scaling about the image centre
    followed by a shift ``(dx, dy)``."""
    height, width = shape
    cx, cy = (width - 1) / 2.0, (height - 1) / 2.0
    a = np.deg2rad(angle_deg)
    c, s = scale * np.cos(a), scale * np.sin(a)
    return np.array(
        [
            [c, -s, cx - c * cx + s * cy + shift[0]],
            [s, c, cy - s * cx - c * cy + shift[1]],
            [0.0, 0.0, 1.0],
        ]
    )


def known_transformations(shape):
    """Ground-truth matrices, one per transformation type, of moderate size."""
    height, width = shape
    affine = rotation_about_center(3.0, shape, scale=1.02, shift=(2.4, -1.7))
    affine[0, 1] += 0.015  # shear
    affine[1, 0] -= 0.01
    bilinear = np.identity(4)
    bilinear[0] = [1.01, 0.012, 4.0e-5, 1.3]
    bilinear[1] = [-0.008, 0.99, -3.0e-5, -0.9]
    return {
        "translation": np.array([[1.0, 0.0, 3.3], [0.0, 1.0, -2.6], [0.0, 0.0, 1.0]]),
        "rigid_body": rotation_about_center(4.0, shape, shift=(2.1, -1.4)),
        "scaled_rotation": rotation_about_center(
            -3.0, shape, scale=1.03, shift=(-1.9, 2.2)
        ),
        "affine": affine,
        "bilinear": bilinear,
    }


def drifting_stack(n_frames=8, shape=(64, 80), kind="rigid_body", seed=3):
    """Stack whose frames drift smoothly, together with the exact matrices.

    Returns ``(stack, tmats)`` with ``stack`` of shape ``(n_frames, *shape)``
    and ``tmats[i]`` mapping frame-0 coordinates onto frame-``i`` coordinates.
    """
    scene = Scene(shape, seed=seed)
    frames, tmats = [], []
    for i in range(n_frames):
        if kind == "translation":
            tmat = np.identity(3)
            tmat[0, 2], tmat[1, 2] = 0.7 * i, -0.45 * i
        else:
            tmat = rotation_about_center(0.6 * i, shape, shift=(0.5 * i, -0.3 * i))
        frames.append(scene.render(tmat))
        tmats.append(tmat)
    return np.stack(frames), np.stack(tmats)
