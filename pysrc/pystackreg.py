# -*- coding: utf-8 -*-
#from pystackreg import stackreg
from . import turboreg
import numpy as np
from tqdm import tqdm

def simple_slice(arr, inds, axis):
    # this does the same as np.take() except only supports simple slicing, not
    # advanced indexing, and thus is much faster
    sl = [slice(None)] * arr.ndim
    sl[axis] = inds
    return arr[sl]

def running_mean(x, N, axis=0):
    pad_width = [[0,0]] * len(x.shape)
    pad_width[axis] = [int(np.ceil(N/2)), int(np.floor(N/2))]
    cumsum = np.cumsum(np.pad(x, pad_width, 'edge'), axis=axis)
    return (cumsum[N:] - cumsum[:-N]) / float(N)

class StackReg:

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
        if transformation not in self._valid_transformations:
            raise Exception("Invalid transformation")
        
        self._transformation = transformation
        self._m = None
        self._tmats = None
        self._refpts = None
        self._movpts = None
    
    def isRegistered(self):
        return self._is_registered

    def register(self, ref, mov):
        self._is_registered = True

        self._m, self._refpts, self._movpts = turboreg._register(ref[:-1, :-1], mov[:-1, :-1], self._transformation)
        return self.get_matrix()
    
    def transform(self, mov, tmat = None):
        if not self.isRegistered():
            raise Exception("Register first")

        if tmat is not None:
            tmat = self._matrix_long_to_short(tmat)
        else:
            tmat = self._m

        return turboreg._transform(mov, tmat)

    def register_transform(self, ref, mov):
        self.register(ref, mov)
        return self.transform(mov)
    
    def get_matrix(self):
        return self._matrix_short_to_long(self._m)

    def set_matrix(self, mat):
        self._m = self._matrix_long_to_short(mat)

    def _matrix_short_to_long(self, m):
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
        return self._refpts, self._movpts

    def transform_stack(self, img, axis=0, tmats=None):
        if tmats is None:
            tmats = self._tmats

        if tmats.shape[0] != img.shape[axis]:
            raise Exception("Number of saved transformation matrices does not match stack length")

        out = img.copy().astype(np.float)

        for i in range(img.shape[axis]):
            slc = [slice(None)] * len(out.shape)
            slc[axis] = slice(i, i+1)
            out[slc] = self.transform(simple_slice(img, i, axis), tmats[i, :, :])

        return out

    def register_stack(self, img, reference='previous', n_frames=1, axis=0, moving_average=1):
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

    def register_transform_stack(self, img, reference='previous', n_frames=1, axis=0, moving_average=1):
        self.register_stack(img, reference, n_frames, axis, moving_average)
        return self.transform_stack(img, axis)

