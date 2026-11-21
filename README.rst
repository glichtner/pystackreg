pyStackReg
==========

.. start-badges

.. image:: https://github.com/glichtner/pystackreg/actions/workflows/ci.yml/badge.svg
    :target: https://github.com/glichtner/pystackreg/actions/workflows/ci.yml
    :alt: Tests

.. image:: https://readthedocs.org/projects/pystackreg/badge/?version=latest
    :target: https://pystackreg.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://badge.fury.io/py/pystackreg.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/pystackreg

.. image:: https://img.shields.io/conda/vn/conda-forge/pystackreg.svg
    :alt: conda-forge package
    :target: https://anaconda.org/conda-forge/pystackreg

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
pyStackReg is used to align (register) one or more images to a common reference image, as is required usually in time-resolved fluorescence or wide-field microscopy. It is directly ported from the source code of the ImageJ plugin ``TurboReg`` and provides additionally the functionality of the ImageJ plugin ``StackReg``, both of which were written by Philippe Thevenaz/EPFL (available at https://bigwww.epfl.ch/thevenaz/turboreg/).

The registration is intensity based: a multi-resolution (pyramid) Marquardt-Levenberg optimizer minimizes the mean squared difference between the reference image and the transformed image, using cubic B-spline interpolation. It works on plain NumPy arrays, needs neither ImageJ nor a Java runtime, and returns the transformation matrices so that they can be stored, inspected and applied to other data.

pyStackReg provides the following five types of distortion:

- translation
- rigid body (translation + rotation)
- scaled rotation (translation + rotation + scaling)
- affine (translation + rotation + scaling + shearing)
- bilinear (non-linear transformation; does not preserve straight lines)

pyStackReg supports the full functionality of StackReg plus some additional options, e.g., using different reference images and having access to the actual transformation matrices (please see the examples below). Note that pyStackReg uses the high quality (i.e. high accuracy) mode of TurboReg that uses cubic spline interpolation for transformation.

Beyond the ImageJ plugins, pyStackReg offers

- three reference strategies for stacks: the previous frame (what StackReg does), the first frame(s), or the mean of the stack, optionally after a moving average,
- multi-threaded registration and transformation of stacks (``n_jobs``),
- memory-efficient processing of large stacks (preallocated, memory-mapped or in-place output),
- transformation of point coordinates with the estimated transformations.

Results agree with the original plugins: for the scaled-rotation, affine and bilinear models the estimated landmarks and transformed images are bit-identical to those of TurboReg in the test suite, for the translation and rigid-body models they agree to better than 1e-6 pixel (see the `validation notes <https://pystackreg.readthedocs.io/en/latest/validation.html>`_).

Please note: The bilinear transformation cannot be propagated, as a combination of bilinear transformations does not generally result in a bilinear transformation. Therefore, stack registration/transform functions won't work with bilinear transformation when using "previous" image as reference image. You can either use another reference ("first" or "mean" for first or mean image, respectively), or try to register/transform each image of the stack separately to its respective previous image (and use the already transformed previous image as reference for the next image).

Scope and known limitations
...........................
- pyStackReg registers 2D images. A 3D array is treated as a *stack of 2D frames* (a time series or a series of sections), each of which gets its own 2D transformation. It does not perform volumetric (3D) registration.
- The optimizer is local. pystackreg (and StackReg/TurboReg) cannot recover large displacements, e.g. rotations by about 30 degrees or more or shifts of a large fraction of the image size (see https://github.com/glichtner/pystackreg/issues/30 for a workaround).
- Both images must show the same kind of contrast; the method is not suited for multi-modal registration.
- Transformed images are interpolated with cubic splines and can contain values slightly outside the range of the input, including negative values. This is the behavior of TurboReg; use ``pystackreg.util.to_uint16`` or an integer ``out`` array to clip and round.
- Images must not contain NaN or infinite values.

Installation
------------
The package is available on conda forge and on PyPi. Binary wheels are provided for Linux (x86-64, aarch64), macOS (x86-64, arm64) and Windows (x86-64) for Python 3.10 and newer; NumPy 1.22 or newer (including NumPy 2) is required.

- Install using **conda**

.. code-block:: console

    conda install pystackreg -c conda-forge

- Install using **pip**

.. code-block:: console

    pip install pystackreg

Installing from source requires a C++11 compiler (``pip install .`` in a checkout of the repository).

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
    # img0.shape: frames x height x width (3D)

    sr = StackReg(StackReg.RIGID_BODY)

    # register 2nd image to 1st
    sr.register(img0[0, :, :], img0[1,:,:])

    # use the transformation from the above registration to register another frame
    out = sr.transform(img1[1,:,:])

The next examples shows how to register and transform a whole stack:

.. code-block:: python

    from pystackreg import StackReg
    from skimage import io

    img0 = io.imread('some_multiframe_image.tif') # 3 dimensions : frames x height x width

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

    img0 = io.imread('some_multiframe_image.tif') # 3 dimensions : frames x height x width
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

    out = np.zeros(img0.shape, dtype=float)

    for i in range(tmats.shape[0]):
        out[i, :, :] = tf.warp(img1[i, :, :], tmats[i, :, :], order=3)

Large stacks can be processed with several threads and without allocating a second, larger copy of the data:

.. code-block:: python

    import numpy as np
    from pystackreg import StackReg

    # memory-mapped 16-bit stack on disk: frames x height x width
    img = np.load('large_stack.npy', mmap_mode='r+')

    sr = StackReg(StackReg.RIGID_BODY)

    # use all CPU cores; results are identical to the sequential run
    tmats = sr.register_stack(img, reference='previous', n_jobs=-1, verbose=True)

    # write the registered frames back into the same array (rounded and
    # clipped to 16 bit) instead of creating a float64 stack of four times the size
    sr.transform_stack(img, n_jobs=-1, out=img)

The transformation matrices map coordinates ``(x, y)`` of the reference image onto coordinates of the unregistered image. To move point coordinates in the same way as the image content is moved, use ``transform_points``:

.. code-block:: python

    sr = StackReg(StackReg.AFFINE)
    tmat = sr.register(ref, mov)

    # positions (x, y) = (column, row) of features detected in the unregistered image
    points = np.array([[120.5, 80.0], [33.2, 14.9]])

    # ... and where they are in the registered image
    registered_points = sr.transform_points(points)

Contributing and support
------------------------
Questions, bug reports and feature requests are welcome in the
`issue tracker <https://github.com/glichtner/pystackreg/issues>`_. For a bug
report, please include the versions of pyStackReg, NumPy and Python and a
minimal example that reproduces the problem.

Pull requests are welcome as well. A development installation compiles the
C++ extension and needs a C++11 compiler::

    git clone https://github.com/glichtner/pystackreg
    cd pystackreg
    pip install -e ".[dev]"
    pytest

Changes to the native code must reproduce the results of earlier releases
(see `Validation <https://pystackreg.readthedocs.io/en/latest/validation.html>`_).
Contributions are accepted under the license of the project
(GPL-3.0-or-later).

Attribution and citation
------------------------
pyStackReg is a port of existing software, not a new registration algorithm:

- **Method:** P. Thévenaz, U.E. Ruttimann, M. Unser: *A Pyramid Approach to
  Subpixel Registration Based on Intensity*, IEEE Transactions on Image
  Processing, vol. 7, no. 1, pp. 27-41, January 1998,
  `doi:10.1109/83.650848 <https://doi.org/10.1109/83.650848>`_.
- **Original software:** the ImageJ plugins
  `TurboReg <https://bigwww.epfl.ch/thevenaz/turboreg/>`_ and
  `StackReg <https://bigwww.epfl.ch/thevenaz/stackreg/>`_ by Philippe Thévenaz,
  Biomedical Imaging Group, École polytechnique fédérale de Lausanne (EPFL).
- **pyStackReg** (translation of TurboReg from Java to C++, Python/NumPy
  interface, stack registration with several reference modes, parallel
  processing): Gregor Lichtner and contributors (see ``AUTHORS.rst``).

If you present or publish results obtained with pyStackReg, please cite the
paper by Thévenaz et al. (1998) - the authors of TurboReg ask for this - and
pyStackReg itself (see ``CITATION.cff``).

License
-------
pyStackReg is free software, distributed under the terms of the
`GNU General Public License <https://www.gnu.org/licenses/gpl-3.0.html>`_,
version 3 or (at your option) any later version; see ``LICENSE``. This is the
license under which TurboReg and StackReg are distributed by their author.

Copyright notices, the conditions of use of the original software and the
license notices of earlier pyStackReg versions (up to 0.2.8, which were
published under different terms) are in ``NOTICE``.
