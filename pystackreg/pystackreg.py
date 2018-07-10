# -*- coding: utf-8 -*-
from . import turboreg
import numpy as np
from tqdm import tqdm

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
    return arr[sl]

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
    pad_width = [[0,0]] * len(x.shape)
    pad_width[axis] = [int(np.ceil(N/2)), int(np.floor(N/2))]
    cumsum = np.cumsum(np.pad(x, pad_width, 'edge'), axis=axis)
    return (cumsum[N:] - cumsum[:-N]) / float(N)

class StackReg:
    """
    Python implementation of the ImageJ/Fiji StackReg plugin ( http://bigwww.epfl.ch/thevenaz/stackreg/ )
    """

    # Transformation
    TRANSLATION = 2
    RIGID_BODY = 3
    SCALED_ROTATION = 4
    AFFINE = 6
    BILINEAR = 8


    _valid_transformations = [
        TRANSLATION, RIGID_BODY, SCALED_ROTATION, AFFINE, BILINEAR
    ]

    _is_registered = False

    def __init__(self, transformation):
        """
        Constructor
        :param transformation: TRANSLATION, RIGID_BODY, SCALED_ROTATION, AFFINE, BILINEAR
        """
        if transformation not in self._valid_transformations:
            raise Exception("Invalid transformation")
        
        self._transformation = transformation
        self._m = None
        self._tmats = None
        self._refpts = None
        self._movpts = None
    
    def is_registered(self):
        """
        Indicates whether register() was already called and a transformation matrix was calculated

        :rtype:  bool
        :return: True if a transformation matrix was already calculated
        """
        return self._is_registered

    def register(self, ref, mov):
        """
        Registers an image to a reference image: Only the transformation matrix
        will be calculated, the image will not be transformed
        (use transform() or register_transform() ).

        :type ref: array_like (Ni..., Nj...)
        :param ref: Reference image (static)

        :type mov: array_like (Ni..., Nj...)
        :param mov:
            The image that should be aligned to the reference image

        :rtype:  ndarray(3,3)
        :return: 2D transformation matrix
        """
        self._is_registered = True

        self._m, self._refpts, self._movpts = turboreg._register(ref[:-1, :-1], mov[:-1, :-1], self._transformation)
        return self.get_matrix()
    
    def transform(self, mov, tmat=None):
        """
        Transform an image according to a previous registration.

        Either a transformation matrix has to be explicitly supplied or
        register() has to be called before calling transform().

        :type mov: array_like (Ni..., Nj...)
        :param mov: The image that will be transformed

        :type tmat:  ndarray(3,3), optional
        :param tmat: The transformation matrix

        :rtype:  ndarray (Ni..., Nj...)
        :return: Transformed image - will be of the same shape as the input image (cropping may occur)
        """
        if tmat is None and not self.is_registered():
            raise Exception("Register first")

        if tmat is not None:
            tmat = self._matrix_long_to_short(tmat)
        else:
            tmat = self._m

        return turboreg._transform(mov, tmat)

    def register_transform(self, ref, mov):
        """
        Register and transform an image to a reference image.

        :type ref: ref: array_like (Ni..., Nj...)
        :param ref: Reference image (static)

        :type mov: array_like (Ni..., Nj...)
        :param mov: The image that should be aligned to the reference image

        :rtype:  ndarray (Ni..., Nj...)
        :return: Transformed image - will be of the same shape as the input image (cropping may occur)
        """
        self.register(ref, mov)
        return self.transform(mov)
    
    def get_matrix(self):
        """
        Get the current transformation matrix

        :rtype:  ndarray(3,3)
        :return: The transformation matrix
        """
        return self._matrix_short_to_long(self._m)

    def set_matrix(self, mat):
        """
        Sets the current transformation matrix

        :type mat: ndarray(3,3)
        :param mat: The transformation matrix
        """
        self._m = self._matrix_long_to_short(mat)

    def _matrix_short_to_long(self, m):
        """
        Converts the transformation matrix from the short form used by TurboReg to
        the canonical form from linear algebra.

        :type m: array_like
        :param m: TurboReg transformation matrix

        :rtype:  ndarray(3,3)
        :return: Canonical transformation matrix
        """
        mat = np.identity(3).astype(np.double)

        if self._transformation == self.TRANSLATION:
            mat[0:2, 2] = m[:,0]
        elif self._transformation in [self.RIGID_BODY, self.SCALED_ROTATION, self.AFFINE]:
            mat[0:2, :] = m[:, [1, 2, 0]]
        elif self._transformation == self.BILINEAR:
            raise Exception("Bilinear transformation matrix not supported")
        else:
            raise Exception("Unexpected transformation")

        return mat

    def _matrix_long_to_short(self, mat):
        """
        Converts the transformation matrix from the canonical form from linear algebra to
        the short form used by TurboReg.

        :type mat: ndarray(3,3)
        :param mat: Canonical transformation matrix

        :rtype:  array_like
        :return: TurboReg transformation matrix
        """
        if self._transformation == self.TRANSLATION:
            m = mat[0:2, 2].reshape((2,1))
        elif self._transformation in [self.RIGID_BODY, self.SCALED_ROTATION, self.AFFINE]:
            m = mat[0:2, [2, 0, 1]]
        elif self._transformation == self.BILINEAR:
            raise Exception("Bilinear transformation matrix not supported")
        else:
            raise Exception("Unexpected transformation")

        return m.astype(np.double)

    def get_points(self):
        """
        Returns the pairs of corresponding points from which the
        transformation matrix can be calculated.

        :rtype:  (ndarray, ndarray)
        :return: (Point coordinates of reference image, Point coordinates of image to be aligned)
        """
        return self._refpts, self._movpts

    def register_stack(self, img, reference='previous', n_frames=1, axis=0, moving_average=1):
        """
        Register a stack of images (movie).
        Note that this function will not transform the image but only calculate the transformation matrices.
        For tranformation use transform_stack() after this function or use register_transform_stack() for a single call.

        :type img: array_like(Ni..., Nj..., Nk...)
        :param img: The 3D stack of images that should be aligned

        :type reference: string, one of ['previous', 'first', 'mean']
        :param reference:
            *  *'previous'*: Aligns each frame (image) to its previous frame in the stack
            *  *'first:'* Aligns each frame (image) to the first frame in the stack - if n_frames is > 1, then
               each frame is aligned to the mean of the first n_frames of the stack
            *  *'mean'*: Aligns each frame (image) to the average of all images in the stack

        :type n_frames: int, optional
        :param n_frames: If reference is 'first', then this parameter specifies the number of frames from the
            beginning of the stack that should be averaged to yield the reference image.

        :type axis: int, optional
        :param axis: The axis of the time dimension

        :type moving_average: int, optional
        :param moving_average:
            If moving_average is greater than 1, a moving average of the stack is first created (using
            a subset size of moving_average) before registration

        :rtype:  ndarray(img.shape[axis], 3, 3)
        :return: The transformation matrix for each image in the stack
        """

        if len(img.shape) != 3:
            raise Exception('Stack must have three dimensions')

        idx_start = 1

        if moving_average > 1:
            idx_start = 0
            size = [0] * len(img.shape)
            size[axis] = moving_average
            img = running_mean(img, moving_average, axis=axis)

        self._tmats = np.repeat(np.identity(3).reshape((1,3,3)),img.shape[axis], axis=0).astype(np.double)

        if reference == 'first':
            ref = np.mean(img.take(range(n_frames), axis=axis), axis=axis)
        elif reference == 'mean':
            ref = img.mean(axis=0)
            idx_start = 0

        for i in tqdm(range(idx_start, img.shape[axis])):

            slc = [slice(None)] * len(img.shape)
            slc[axis] = slice(i, i+1)

            if reference == 'previous':
                ref = img.take(i-1, axis=axis)

            self._tmats[i, :, :] = self.register(ref, simple_slice(img, i, axis))

            if reference == 'previous' and i > 0:
                self._tmats[i, :, :] = np.matmul(self._tmats[i, :, :], self._tmats[i-1, :, :])

        return self._tmats

    def transform_stack(self, img, axis=0, tmats=None):
        """
        Transform a stack after registration.

        :type img: array_like(Ni..., Nj..., Nk...)
        :param img: The 3D stack of images that should be aligned

        :type axis: int, optional
        :param axis: The axis of the time dimension

        :type tmats: ndarray(img.shape[axis], 3, 3), optional
        :param tmats: The transformation matrix for each image in the stack

        :rtype:  ndarray(Ni..., Nj..., Nk...)
        :return: The transformed stack
        """
        if tmats is None:
            tmats = self._tmats

        if tmats.shape[0] != img.shape[axis]:
            raise Exception("Number of saved transformation matrices does not match stack length")

        out = img.copy().astype(np.float)

        for i in range(img.shape[axis]):
            slc = [slice(None)] * len(out.shape)
            slc[axis] = slice(i, i + 1)
            out[slc] = self.transform(simple_slice(img, i, axis), tmats[i, :, :])

        return out

    def register_transform_stack(self, img, reference='previous', n_frames=1, axis=0, moving_average=1):
        """
        Register and transform stack of images (movie).

        :type img: array_like(Ni..., Nj..., Nk...)
        :param img: The 3D stack of images that should be aligned

        :type reference: string, one of ['previous', 'first', 'mean']
        :param reference:
            *  *'previous'*: Aligns each frame (image) to its previous frame in the stack
            *  *'first:'* Aligns each frame (image) to the first frame in the stack - if n_frames is > 1, then
               each frame is aligned to the mean of the first n_frames of the stack
            *  *'mean'*: Aligns each frame (image) to the average of all images in the stack

        :type n_frames: int, optional
        :param n_frames: If reference is 'first', then this parameter specifies the number of frames from the
            beginning of the stack that should be averaged to yield the reference image.

        :type axis: int, optional
        :param axis: The axis of the time dimension

        :type moving_average: int, optional
        :param moving_average: If moving_average is greater than 1, a moving average of the stack is first created (using
            a subset size of moving_average) before registration

        :rtype:  ndarray(Ni..., Nj..., Nk...)
        :return: The transformed stack
        """
        self.register_stack(img, reference, n_frames, axis, moving_average)
        return self.transform_stack(img, axis)

