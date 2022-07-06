"""
Collection of utility functions for pystackreg
"""
import numpy as np


def to_uint16(arr):
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


def to_int_dtype(arr: np.ndarray, dtype: type):
    """
    Converts image data array to an integer array of the given dtype.

    .. warning:: Data loss may occur as data points outside of the datatypes min/max are
       clipped.

    StackReg generally output float arrays. This utility function can be used to convert
    such an output array to integer data by clipping to the bounds of the integer dtype
    and rounding the floats.

    :type arr: array_like (Ni..., M, Nk...)
    :param arr: Source array (usually float)

    :type dtype: type
    :param dtype: Integer datatype (e.g. np.uint16)

    :rtype:  ndarray(Ni..., Nj..., Nk...)
    :return: Input array clipped & rounded to integer dtype
    """
    assert type(arr) == np.ndarray, "Input must be a numpy array"

    if not np.issubdtype(dtype, np.integer):
        return arr

    ii = np.iinfo(dtype)

    return np.round(arr.clip(min=ii.min, max=ii.max)).astype(dtype)


def simple_slice(arr, inds, axis):
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

    sl = [slice(None)] * arr.ndim
    sl[axis] = inds
    return arr[tuple(sl)]


def running_mean(x, N, axis=0):
    """
    Calculate running mean (=moving average) across a given axis.

    The array is padded with the first and last value such that
    the resulting running mean has the same dimensions as the input array.

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
    pad_width = [[0, 0]] * len(x.shape)
    pad_width[axis] = [int(np.ceil(N / 2)), int(np.floor(N / 2))]
    cumsum = np.cumsum(np.pad(x, pad_width, "edge"), axis=axis)
    return (cumsum[N:] - cumsum[:-N]) / float(N)
