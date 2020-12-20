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

    assert type(arr) == np.ndarray, "Input must be a numpy array"

    return np.round(arr.clip(min=0, max=65535)).astype(np.uint16)
