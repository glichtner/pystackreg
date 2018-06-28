==================================================
pystackreg: A python/c++ port of TurboReg/Stackreg
==================================================

|Pypi Package| |Build Status| |Coverage Status| |Code Health| |License|

Summary
-------
Python/C++ port of the ImageJ extension TurboReg/StackReg written by Philippe Thevenaz/EPFL



Installation
------------
The package is available on PyPi. Install it using:

.. code-block:: python

    pip install pystackreg



Usage
-----
The following example opens two different files and registers them using all different possible transformations

.. code-block:: python

    from pystackreg import StackReg
    from skimage import io
    
    ref = io.imread('some_original_image.tif')
    mov = io.imread('some_changed_image.tif')
    
    sr = StackReg(StackReg.TRANSLATION)
    out_tra = sr.register_transform(ref, mov)
    
    sr = StackReg(StackReg.RIGID_BODY)
    out_rot = sr.register_transform(ref, mov)
    
    sr = StackReg(StackReg.SCALED_ROTATION)
    out_sca = sr.register_transform(ref, mov)
    
    sr = StackReg(StackReg.AFFINE)
    out_aff = sr.register_transform(ref, mov)



The next example shows how to separate registration from transformation (e.g., to register in one color channel and then use that information to transform another color channel):


.. code-block:: python

    from pystackreg import StackReg
    from skimage import io
    
    img0 = io.imread('some_multiframe_image.tif') # 3 dimensions : frames x width x height
    img1 = io.imread('another_multiframe_image.tif')
    
    sr = StackReg(StackReg.RIGID_BODY)
    sr.register(img0[0, :, :], img0[1,:,:]) # register 2nd image to 1st
    out = sr.transform(img1[1,:,:]) # use the transformation from the above registration to register another frame

The next examples shows how to register and transform a whole stack:

.. code-block:: python

    from pystackreg import StackReg
    from skimage import io
    
    img0 = io.imread('some_multiframe_image.tif') # 3 dimensions : frames x width x height
    
    sr = StackReg(StackReg.RIGID_BODY)
    
    # register to first image
    out_first = register_transform_stack(img0, reference='first')
    
    # register to mean image
    out_mean = register_transform_stack(img0, reference='mean')
    
    # register each frame to the previous (already registered) one 
    out_previous = register_transform_stack(img0, reference='previous')
    
    # register to mean of first 10 images
    out_first10 = register_transform_stack(img0, reference='first', n_frames=10)
    
    # calculate a moving average of 10 images, then register the moving average to the mean of the first 10 images
    # and transform the original image (not the moving average)
    out_moving10 = register_transform_stack(img0, reference='first', n_frames=10, moving_average = 10)


Author information
-------------------
This is a port of the original Java code by Philippe Thevenaz to C++ with a Python wrapper around it. All credit goes to the original author:
::

    /*====================================================================
    | Philippe Thevenaz
    | EPFL/STI/IMT/LIB/BM.4.137
    | Station 17
    | CH-1015 Lausanne VD
    | Switzerland
    |
    | phone (CET): +41(21)693.51.61
    | fax: +41(21)693.37.01
    | RFC-822: philippe.thevenaz@epfl.ch
    | X-400: /C=ch/A=400net/P=switch/O=epfl/S=thevenaz/G=philippe/
    | URL: http://bigwww.epfl.ch/
    \===================================================================*/
    
    /*====================================================================
    | This work is based on the following paper:
    |
    | P. Thevenaz, U.E. Ruttimann, M. Unser
    | A Pyramid Approach to Subpixel Registration Based on Intensity
    | IEEE Transactions on Image Processing
    | vol. 7, no. 1, pp. 27-41, January 1998.
    |
    | This paper is available on-line at
    | http://bigwww.epfl.ch/publications/thevenaz9801.html
    |
    | Other relevant on-line publications are available at
    | http://bigwww.epfl.ch/publications/
    \===================================================================*/

License
-------
Below is the license of TurboReg/StackReg:

::

    /*====================================================================
    | Additional help available at http://bigwww.epfl.ch/thevenaz/turboreg/
    |
    | You'll be free to use this software for research purposes, but you
    | should not redistribute it without our consent. In addition, we expect
    | you to include a citation or acknowledgment whenever you present or
    | publish results that are based on it.
    \===================================================================*/
    

