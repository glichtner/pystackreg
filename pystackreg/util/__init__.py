"""
Collection of utility functions for pystackreg
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

__all__ = [
    "max_displacement",
    "running_mean",
    "simple_slice",
    "to_int_dtype",
    "to_uint16",
]


def to_uint16(arr: npt.ArrayLike) -> np.ndarray:
    """
    Converts image data array to a 16 bit unsigned integer array.

    .. warning:: Data loss may occur as data points outside of 16 bit (0-65535) are
       clipped.

    StackReg generally output float arrays. This utility function can be used to convert
    such an output array to 16 bit data by clipping to 16 bit (to remove values < 0 and
    greater 16 bit) and rounding the floats.

    :type arr: array_like (Ni..., M, Nk...)
    :param arr: Source array (usually float)

    :rtype:  ndarray(Ni..., Nj..., Nk...)
    :return: Input array clipped & rounded to unsigned 16 bit integer
    """
    return to_int_dtype(arr, np.uint16)


def to_int_dtype(arr: npt.ArrayLike, dtype: npt.DTypeLike) -> np.ndarray:
    """
    Converts image data array to an integer array of the given dtype.

    .. warning:: Data loss may occur as data points outside of the datatypes min/max are
       clipped.

    StackReg generally output float arrays. This utility function can be used to convert
    such an output array to integer data by clipping to the bounds of the integer dtype
    and rounding the floats (half to even, as :func:`numpy.round` does).

    :type arr: array_like (Ni..., M, Nk...)
    :param arr: Source array (usually float)

    :type dtype: type
    :param dtype: Integer datatype (e.g. np.uint16). For a non-integer dtype the
                  input is returned unchanged.

    :rtype:  ndarray(Ni..., Nj..., Nk...)
    :return: Input array clipped & rounded to integer dtype
    """
    arr = np.asarray(arr)

    if not np.issubdtype(dtype, np.integer):
        return arr

    ii = np.iinfo(dtype)

    return np.round(arr.clip(min=ii.min, max=ii.max)).astype(dtype)


def simple_slice(arr: np.ndarray, inds: Any, axis: int) -> np.ndarray:
    """
    Take elements from an array along an axis.

    This does the same as np.take() except only supports simple slicing, not
    advanced indexing, and thus is much faster

    :type arr: array_like (Ni..., M, Nk...)
    :param arr: The source array to slice from

    :type inds: int or array_like (Nj...)
    :param inds:
        The indices of the values to extract

    :type axis: int
    :param axis: The axis over which to select values

    :rtype:  ndarray(Ni..., Nj..., Nk...)
    :return: The returned array has the same type as arr
    """

    sl: list[Any] = [slice(None)] * arr.ndim
    sl[axis] = inds
    return arr[tuple(sl)]


def running_mean(x: npt.ArrayLike, N: int, axis: int = 0) -> np.ndarray:
    """
    Calculate running mean (=moving average) across a given axis.

    The array is padded with the first and last value such that
    the resulting running mean has the same dimensions as the input array.

    .. note:: This function builds padded and cumulated copies of the whole
       array. :meth:`pystackreg.StackReg.register_stack` computes the same
       moving average frame by frame and does not use it.

    :type x: array_like (Ni..., Nj..., Nk...)
    :param x: The source array

    :type N: int
    :param N:
        Number of elements to average over

    :type axis: int, optional
    :param axis: The axis across which the running mean is calculated

    :rtype:  ndarray(Ni..., Nj..., Nk...)
    :return: The returned array has the same shape and type as x
    """
    x = np.asarray(x)
    pad_width = [[0, 0]] * len(x.shape)
    pad_width[axis] = [int(np.ceil(N / 2)), int(np.floor(N / 2))]
    cumsum = np.cumsum(np.pad(x, pad_width, "edge"), axis=axis)
    upper = simple_slice(cumsum, slice(N, None), axis)
    lower = simple_slice(cumsum, slice(None, -N), axis)
    return (upper - lower) / float(N)


def max_displacement(
    tmats: npt.ArrayLike, shape: tuple[int, int]
) -> np.ndarray | float:
    """
    Largest distance (in pixels) by which a transformation moves any point of
    an image.

    This is a single, interpretable measure of the *amplitude* of a
    transformation that works for all transformation types. It can be used
    to detect implausible registration results, e.g. frames for which an
    affine registration converged to a strongly scaled or sheared solution::

        tmats = sr.register_stack(stack, reference="first")
        suspicious = max_displacement(tmats, stack.shape[1:]) > 25

    The displacement of an affine or bilinear map over a rectangle is largest
    at one of the corners, so only the four image corners are evaluated.

    :type tmats: ndarray(3, 3), ndarray(4, 4) or a stack of them (N, 3, 3) /
                 (N, 4, 4)
    :param tmats: Transformation matrix or matrices as returned by
                  :meth:`pystackreg.StackReg.register` and
                  :meth:`pystackreg.StackReg.register_stack`

    :type shape: (int, int)
    :param shape: Image shape ``(rows, columns)``

    :rtype:  float or ndarray(N)
    :return: Maximal displacement in pixels (one value per matrix)
    """
    mats = np.asarray(tmats, dtype=np.float64)
    single = mats.ndim == 2
    if single:
        mats = mats[np.newaxis]
    if mats.ndim != 3 or mats.shape[1:] not in ((3, 3), (4, 4)):
        raise ValueError(
            "tmats must have shape (3, 3), (4, 4), (N, 3, 3) or (N, 4, 4), "
            f"got {np.shape(tmats)}"
        )
    rows, cols = shape
    x = np.array([0.0, cols - 1.0, 0.0, cols - 1.0])
    y = np.array([0.0, 0.0, rows - 1.0, rows - 1.0])
    if mats.shape[1] == 3:
        vec = np.stack([x, y, np.ones(4)])
    else:
        vec = np.stack([x, y, x * y, np.ones(4)])
    mapped = mats[:, :2, :] @ vec  # (N, 2, 4)
    dist = np.hypot(mapped[:, 0, :] - x, mapped[:, 1, :] - y).max(axis=1)
    return float(dist[0]) if single else dist
