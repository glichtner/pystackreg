# -*- coding: utf-8 -*-
#from pystackreg import stackreg
from . import stackreg

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
    

