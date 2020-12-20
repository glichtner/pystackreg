pyStackReg
==========

.. start-badges

.. image:: https://ci.appveyor.com/api/projects/status/3kqq8qyc9b7o1coe?svg=true
    :target: https://ci.appveyor.com/api/projects/status/3kqq8qyc9b7o1coe?svg=true
    :alt: Build status

.. image:: https://readthedocs.org/projects/pystackreg/badge/?version=latest
    :target: https://pystackreg.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://badge.fury.io/py/pystackreg.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/pystackreg

.. image:: https://img.shields.io/pypi/pyversions/pystackreg.svg
    :alt: Supported Python Versions
    :target: https://pypi.org/project/pystackreg/

.. image:: https://pepy.tech/badge/pystackreg
    :alt: Downloads
    :target: https://pepy.tech/project/pystackreg/

.. end-badges





Summary
-------
Python/C++ port of the ImageJ extension TurboReg/StackReg written by Philippe Thevenaz/EPFL.

A python extension for the automatic alignment of a source image or a stack (movie) to a target image/reference frame.

Description
-----------
pyStackReg is used to align (register) one or more images to a common reference image, as is required usually in time-resolved fluorescence or wide-field microscopy. It is directly ported from the source code of the ImageJ plugin ``TurboReg`` and provides additionally the functionality of the ImageJ plugin ``StackReg``, both of which were written by Philippe Thevenaz/EPFL (available at http://bigwww.epfl.ch/thevenaz/turboreg/).

pyStackReg provides the following five types of distortion:

- translation
- rigid body (translation + rotation)
- scaled rotation (translation + rotation + scaling)
- affine (translation + rotation + scaling + shearing)
- bilinear (non-linear transformation; does not preserve straight lines)

pyStackReg supports the full functionality of StackReg plus some additional options, e.g., using different reference images and having access to the actual transformation matrices (please see the examples below).

Please note: The bilinear transformation cannot be propagated, as a combination of bilinear transformations does not generally result in a bilinear transformation. Therefore, stack registration/transform functions won't work with bilinear transformation when using "previous" image as reference image. You can either use another reference ("first" or "mean" for first or mean image, respectively), or try to register/transform each image of the stack separately to its respective previous image (and use the already transformed previous image as reference for the next image).


Installation
------------
The package is available on conda forge and on PyPi.

- Install using **conda**

.. code-block:: python

    conda install pystackreg -c conda-forge

- Install using **pip**

.. code-block:: python

    pip install pystackreg


Documentation
-------------
The documentation can be found on readthedocs:

https://pystackreg.readthedocs.io/

Tutorial
--------
* A tutorial notebook can be found in the `examples/notebooks` folder
  or statically here: https://pystackreg.readthedocs.io/en/latest/tutorial.html

Usage
-----
The following example opens two different files and registers them using all different possible transformations

.. code-block:: python

    from pystackreg import StackReg
    from skimage import io

    #load reference and "moved" image
    ref = io.imread('some_original_image.tif')
    mov = io.imread('some_changed_image.tif')

    #Translational transformation
    sr = StackReg(StackReg.TRANSLATION)
    out_tra = sr.register_transform(ref, mov)

    #Rigid Body transformation
    sr = StackReg(StackReg.RIGID_BODY)
    out_rot = sr.register_transform(ref, mov)

    #Scaled Rotation transformation
    sr = StackReg(StackReg.SCALED_ROTATION)
    out_sca = sr.register_transform(ref, mov)

    #Affine transformation
    sr = StackReg(StackReg.AFFINE)
    out_aff = sr.register_transform(ref, mov)

    #Bilinear transformation
    sr = StackReg(StackReg.BILINEAR)
    out_bil = sr.register_transform(ref, mov)


The next example shows how to separate registration from transformation (e.g., to register in one color channel and then use that information to transform another color channel):


.. code-block:: python

    from pystackreg import StackReg
    from skimage import io

    img0 = io.imread('some_multiframe_image.tif')
    img1 = io.imread('another_multiframe_image.tif')
    # img0.shape: frames x width x height (3D)

    sr = StackReg(StackReg.RIGID_BODY)

    # register 2nd image to 1st
    sr.register(img0[0, :, :], img0[1,:,:])

    # use the transformation from the above registration to register another frame
    out = sr.transform(img1[1,:,:])

The next examples shows how to register and transform a whole stack:

.. code-block:: python

    from pystackreg import StackReg
    from skimage import io

    img0 = io.imread('some_multiframe_image.tif') # 3 dimensions : frames x width x height

    sr = StackReg(StackReg.RIGID_BODY)

    # register each frame to the previous (already registered) one
    # this is what the original StackReg ImageJ plugin uses
    out_previous = sr.register_transform_stack(img0, reference='previous')

    # register to first image
    out_first = sr.register_transform_stack(img0, reference='first')

    # register to mean image
    out_mean = sr.register_transform_stack(img0, reference='mean')

    # register to mean of first 10 images
    out_first10 = sr.register_transform_stack(img0, reference='first', n_frames=10)

    # calculate a moving average of 10 images, then register the moving average to the mean of
    # the first 10 images and transform the original image (not the moving average)
    out_moving10 = sr.register_transform_stack(img0, reference='first', n_frames=10, moving_average = 10)

The next example shows how to separate registration from transformation for a stack (e.g., to register in one color channel and then use that information to transform another color channel):

.. code-block:: python

    from pystackreg import StackReg
    from skimage import io

    img0 = io.imread('some_multiframe_image.tif') # 3 dimensions : frames x width x height
    img1 = io.imread('another_multiframe_image.tif') # same shape as img0

    # both stacks must have the same shape
    assert img0.shape == img1.shape

    sr = StackReg(StackReg.RIGID_BODY)

    # register each frame to the previous (already registered) one
    # this is what the original StackReg ImageJ plugin uses
    tmats = sr.register_stack(img0, reference='previous')
    out = sr.transform_stack(img1)

    # tmats contains the transformation matrices -> they can be saved
    # and loaded at another time
    import numpy as np
    np.save('transformation_matrices.npy', tmats)

    tmats_loaded = np.load('transformation_matrices.npy')

    # make sure you use the correct transformation here!
    sr = StackReg(StackReg.RIGID_BODY)

    # transform stack using the tmats loaded from file
    sr.transform_stack(img1, tmats=tmats_loaded)

    # with the transformation matrices at hand you can also
    # use the transformation algorithms from other packages:
    from skimage import transform as tf

    out = np.zeros(img0.shape).astype(np.float)

    for i in range(tmats.shape[0]):
        tform = tf.AffineTransform(matrix=tmats[i, :, :])
        out[i, :, :] = tf.warp(img1[i, :, :], tform)


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

::

    You are free to use this software for commercial and non-commercial
    purposes. However, we expect you to include a citation or acknowledgement
    whenever you present or publish research results that are based
    on this software. You are free to modify this software or derive
    works from it, but you are only allowed to distribute it under the
    same terms as this license specifies. Additionally, you must include
    a reference to the research paper above in all software and works
    derived from this software.
