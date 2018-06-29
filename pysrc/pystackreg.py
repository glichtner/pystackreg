# -*- coding: utf-8 -*-
#from pystackreg import stackreg
from . import stackreg
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

    _valid_transformations = [
        TRANSLATION, RIGID_BODY, SCALED_ROTATION, AFFINE
    ]

    _is_registered = False

    def __init__(self, transformation):
        if not transformation in self._valid_transformations:
            raise Exception("Invalid transformation")
        
        self._transformation = transformation
    
    def isRegistered(self):
        return self._is_registered

    def register(self, ref, mov):
        self._is_registered = True
        self.m, self.refpts, self.movpts = stackreg._register(ref, mov, self._transformation)
    
    def transform(self, mov):
        if not self.isRegistered():
            raise Exception("Register first")
        return stackreg._transform(mov, self.m)

    def register_transform(self, ref, mov):
        self.register(ref, mov)
        return self.transform(mov)
    
    def get_matrix(self):
        return self.m

    def get_points(self):
        return (self.refpts, self.movpts)

    def register_transform_stack(self, img, reference='first', n_frames=1, axis=0, moving_average=1):
        if len(img.shape) != 3:
            raise Exception('Stack must have three dimensions')
        
        out = img.copy().astype(np.float)
        
        idx_start = 1
        
        if moving_average > 1:
            from scipy.ndimage import uniform_filter
            idx_start = 0
            size = [0] * len(img.shape)
            size[axis] = moving_average
            img_orig = img.copy()
            img = running_mean(img, moving_average, axis=axis)
            
        if reference == 'first':
            ref = np.mean(img.take(range(n_frames), axis=axis), axis=axis)
        elif reference == 'mean':
            ref = img.mean(axis=0)
            idx_start = 0
            
      
        for i in tqdm(range(idx_start, img.shape[axis])):
            
            slc = [slice(None)] * len(out.shape)
            slc[axis] = slice(i,i+1)
            
            if reference == 'previous':
                ref = out.take(i-1, axis=axis)
            
            if moving_average > 1:
                self.register(ref, simple_slice(img, i, axis))
                out[slc] = self.transform(simple_slice(img_orig, i, axis))
            else:
                out[slc] = self.register_transform(ref, simple_slice(img, i, axis))
            
        
        return out
        
        

