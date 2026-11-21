User guide
==========

This page explains the concepts behind :class:`pystackreg.StackReg`: what the
five transformation types mean, how the returned matrices are to be read, how
stacks are registered, and where the limits of the method are.

What pystackreg does
--------------------

Given a *reference* image and a *moving* image of the same size, pystackreg
estimates the geometric transformation that aligns the moving image with the
reference (:meth:`~pystackreg.StackReg.register`) and resamples the moving
image accordingly (:meth:`~pystackreg.StackReg.transform`).

The estimation is the algorithm of TurboReg [Thevenaz1998]_: the mean squared
intensity difference between the two images is minimized with a
Marquardt-Levenberg optimizer, coarse-to-fine on an image pyramid, with images
modelled by cubic B-splines. pystackreg always uses TurboReg's "accurate" mode.
The result has sub-pixel precision; on well-textured test images the estimated
transformations deviate from the ground truth by less than a hundredth of a
pixel (see :doc:`validation`).

.. note::

   pystackreg registers **2D images**. A 3D array is treated as a *stack of 2D
   frames* - a time series or a series of sections - and every frame gets its
   own 2D transformation. This is different from registering two *volumes*
   with a 3D transformation, which pystackreg does not do.

Transformation types
--------------------

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Constant
     - Parameters
     - Motion model
   * - ``StackReg.TRANSLATION``
     - 2
     - shift
   * - ``StackReg.RIGID_BODY``
     - 3
     - shift and rotation
   * - ``StackReg.SCALED_ROTATION``
     - 4
     - shift, rotation and isotropic scaling
   * - ``StackReg.AFFINE``
     - 6
     - shift, rotation, anisotropic scaling and shear; straight lines stay
       straight and parallel lines stay parallel
   * - ``StackReg.BILINEAR``
     - 8
     - non-linear; straight lines are in general not preserved

Choose the model with the fewest parameters that can describe the motion in
your data. Models with more parameters can fit noise and are more likely to
converge to an implausible solution (see :ref:`amplitude`).

.. _matrix-convention:

Coordinates and transformation matrices
---------------------------------------

Coordinates are written :math:`(x, y)` with :math:`x` the **column** index and
:math:`y` the **row** index of a pixel; the origin is the centre of the
top-left pixel. An image ``img`` is therefore sampled at ``img[y, x]``.

:meth:`~pystackreg.StackReg.register` returns a matrix :math:`T` that maps a
position in the **reference** image onto the corresponding position in the
**moving** image:

.. math::

   \begin{pmatrix} x_\mathrm{mov} \\ y_\mathrm{mov} \\ 1 \end{pmatrix}
   = T \begin{pmatrix} x_\mathrm{ref} \\ y_\mathrm{ref} \\ 1 \end{pmatrix}

This is the direction needed for resampling:
:meth:`~pystackreg.StackReg.transform` computes, for every pixel
:math:`(x_\mathrm{ref}, y_\mathrm{ref})` of the output, the value of the
moving image at :math:`T (x_\mathrm{ref}, y_\mathrm{ref}, 1)^\top`.

The matrices of the four linear models are 3 x 3:

.. math::

   T_\mathrm{translation} =
   \begin{pmatrix} 1 & 0 & t_x \\ 0 & 1 & t_y \\ 0 & 0 & 1 \end{pmatrix},
   \quad
   T_\mathrm{rigid} =
   \begin{pmatrix} \cos\theta & -\sin\theta & t_x \\
                   \sin\theta & \cos\theta & t_y \\ 0 & 0 & 1 \end{pmatrix},
   \quad
   T_\mathrm{scaled} =
   \begin{pmatrix} a & -b & t_x \\ b & a & t_y \\ 0 & 0 & 1 \end{pmatrix},
   \quad
   T_\mathrm{affine} =
   \begin{pmatrix} a_{11} & a_{12} & t_x \\ a_{21} & a_{22} & t_y \\
                   0 & 0 & 1 \end{pmatrix}

For the scaled rotation the scale factor is :math:`\sqrt{a^2 + b^2}` and the
angle :math:`\operatorname{atan2}(b, a)`.

.. note::

   Rotation and scaling act about the **origin** (the top-left pixel), not
   about the image centre. :math:`(t_x, t_y)` is therefore *not* the shift of
   the image content when the matrix contains a rotation. To obtain the
   displacement of a particular point, e.g. the image centre :math:`c`,
   compute :math:`T c - c`::

       rows, cols = img.shape
       centre = np.array([(cols - 1) / 2, (rows - 1) / 2, 1.0])
       shift_of_centre = (tmat @ centre)[:2] - centre[:2]
       angle = np.degrees(np.arctan2(tmat[1, 0], tmat[0, 0]))

The **bilinear** model is

.. math::

   x_\mathrm{mov} = a_1 x + a_2 y + a_3 x y + a_0, \qquad
   y_\mathrm{mov} = b_1 x + b_2 y + b_3 x y + b_0

and is returned as a 4 x 4 matrix acting on :math:`(x, y, xy, 1)^\top`; its
first two rows are :math:`(a_1, a_2, a_3, a_0)` and
:math:`(b_1, b_2, b_3, b_0)`, the other two rows only pad the matrix. Bilinear
transformations cannot be chained by matrix multiplication.

Landmarks
.........

TurboReg parametrizes a transformation by pairs of corresponding points
("landmarks") in the two images; the matrix is derived from them.
:meth:`~pystackreg.StackReg.get_points` returns the landmarks of the reference
image (fixed positions that depend only on the image size) and the refined
landmarks of the moving image. There are 1, 3, 2, 3 and 4 pairs for the
translation, rigid-body, scaled-rotation, affine and bilinear model.

Transforming point coordinates
..............................

To move point coordinates - object centroids, contours, regions of interest -
in the same way as :meth:`~pystackreg.StackReg.transform` moves the image
content, use :meth:`~pystackreg.StackReg.transform_points`::

    sr = StackReg(StackReg.AFFINE)
    tmat = sr.register(ref, mov)

    points = np.array([[120.5, 80.0], [33.2, 14.9]])   # (x, y) in `mov`
    registered = sr.transform_points(points)           # (x, y) after registration

This applies the *inverse* of the matrix, because the matrix points from the
reference to the moving image. Multiplying coordinates by the matrix itself
(``inverse=True``) answers the opposite question: where in the moving image is
a given position of the reference? Coordinates obtained from NumPy indexing
(``np.argwhere``, ``skimage.measure.regionprops``) are in (row, column) order
and have to be flipped: ``sr.transform_points(rc[:, ::-1])[:, ::-1]``.

For a stack, apply the matrix of each frame to the points of that frame::

    tmats = sr.register_stack(stack, reference="first")
    registered = [sr.transform_points(p, t) for p, t in zip(points_per_frame, tmats)]

Using the matrices with other libraries
.......................................

The 3 x 3 matrices can be applied with other libraries, for example to use a
different interpolation. With SciPy, which works in (row, column) order::

    from scipy import ndimage
    rc = tmat[[1, 0, 2]][:, [1, 0, 2]]            # swap x and y
    out = ndimage.affine_transform(mov, rc[:2, :2], offset=rc[:2, 2], order=3)

``order=3`` is the same cubic B-spline interpolation that pystackreg uses;
away from the border the result equals that of
:meth:`~pystackreg.StackReg.transform` to single precision. (The two differ in
how they treat positions outside the image: pystackreg writes 0 as soon as the
nearest source pixel lies outside.)

The matrices also follow the inverse-map convention of scikit-image and can be
passed to ``skimage.transform.warp`` directly::

    from skimage import transform
    out = transform.warp(mov, tmat, order=1, preserve_range=True)

Note that scikit-image's ``order=3`` is a different cubic interpolation scheme
than the B-spline interpolation of TurboReg; on a textured test image the two
differ by about 1 % of the intensity range. ``order=1`` (bilinear
interpolation) never leaves the range of the input values, at the price of
some smoothing.

Registering stacks
------------------

:meth:`~pystackreg.StackReg.register_stack` computes one matrix per frame,
:meth:`~pystackreg.StackReg.transform_stack` applies them, and
:meth:`~pystackreg.StackReg.register_transform_stack` does both. The frame
axis is the first axis by default (``axis=0``, i.e. the stack has the shape
frames x rows x columns). pystackreg warns if another axis looks more like the
frame axis; pass ``suppress_axis_warning=True`` if you are sure.

Reference strategies
....................

``reference="previous"`` (default; what the StackReg plugin does)
    Every frame is registered to the frame before it and the pairwise
    transformations are chained: :math:`T_i = P_i\,T_{i-1}` with
    :math:`T_0 = I`. All frames end up in the coordinate system of the first
    frame. Consecutive frames are most similar, so this is robust for slowly
    changing image content. Small errors of the pairwise registrations
    accumulate along the stack. Not available for the bilinear model.

``reference="first"``
    Every frame is registered directly to the first frame, or to the mean of
    the first ``n_frames`` frames. No error accumulation, but the displacement
    of late frames relative to the first frame must stay within the capture
    range of the method, and the image content must not change much.

``reference="mean"``
    Every frame is registered to the mean of all frames. The mean image has
    the best signal-to-noise ratio but is blurred if the stack moves a lot;
    it then helps to register twice (transform, then register the result to
    its sharper mean).

``moving_average=N``
    Registers a moving average over ``N`` frames instead of the individual
    frames, which helps with very noisy data. The resulting matrices are
    meant to be applied to the original frames.

The matrices can be stored (``np.save``) and applied later, or to another
stack of the same length, e.g. a second colour channel::

    tmats = sr.register_stack(channel_1, reference="previous")
    registered_2 = sr.transform_stack(channel_2, tmats=tmats)

For large stacks see :doc:`performance`.

.. _amplitude:

Checking the plausibility of a result
-------------------------------------

The optimizer always returns a transformation, also when the images do not
match. Models with many parameters (affine, bilinear) can then produce strong
scaling or shear. :func:`pystackreg.util.max_displacement` condenses a matrix
into one number, the largest distance in pixels by which it moves any point of
the image, which makes implausible results easy to find::

    from pystackreg.util import max_displacement

    tmats = sr.register_stack(stack, reference="first")
    displacement = max_displacement(tmats, stack.shape[1:])
    suspicious = np.flatnonzero(displacement > 30)      # threshold in pixels

What to do with such frames depends on the data: register them with a more
constrained model (e.g. ``RIGID_BODY``), replace their matrix by that of a
neighbouring frame, or exclude them. pystackreg does not limit the
transformation during the optimization, because this would change the
behaviour of the TurboReg algorithm.

.. _limitations:

Limitations
-----------

Limited capture range
    The optimizer is local. Starting from the identity, it finds the nearest
    minimum of the image difference. On the test scenes rotations up to about
    20 degrees and shifts up to about a quarter of the image size were
    recovered, 30 degrees and more were not; the limits depend on the image
    content. For larger displacements, pre-align the images coarsely (e.g.
    try rotations by multiples of 90 degrees and keep the best match, or
    estimate a shift by cross-correlation) and let pystackreg refine the
    result.

Same contrast in both images
    The criterion is the squared intensity difference. Images of different
    modalities, with inverted contrast, or with a strongly different
    brightness cannot be registered directly. A linear rescaling of *both*
    images does not change the result.

Structure is required
    Images without structure (constant images) carry no information about
    the transformation; the result is the identity or meaningless.

Interpolation overshoot and negative values
    Transformed images are computed with cubic spline interpolation. Next to
    sharp edges the result can lie slightly outside the range of the input,
    so an image without negative values can acquire some. This is the
    behaviour of TurboReg and is not corrected silently. Convert with
    :func:`pystackreg.util.to_uint16` / :func:`pystackreg.util.to_int_dtype`
    or pass an integer ``out`` array to
    :meth:`~pystackreg.StackReg.transform_stack` to round and clip. If you
    need a different interpolation, apply the matrices with another library
    (see above).

Borders
    Output pixels whose source position lies outside the moving image are set
    to 0. The last row and the last column of the images are not used for the
    registration (as in StackReg).

Non-finite values
    Images with NaN or infinite values are rejected by the registration. In
    :meth:`~pystackreg.StackReg.transform` a single NaN spreads over the
    whole output, because the spline coefficients are computed with a
    recursive filter. Replace such pixels before calling pystackreg.

Image size
    Both images must have the same shape and at least 3 x 3 pixels. Images
    smaller than about 24 pixels along one side are registered at full
    resolution only, without a pyramid, which reduces the capture range.

Colour images
    Register one channel (or a grey-value combination) and apply the
    matrices to all channels.

.. [Thevenaz1998] P. Thévenaz, U.E. Ruttimann, M. Unser, "A Pyramid Approach
   to Subpixel Registration Based on Intensity", IEEE Transactions on Image
   Processing 7(1), 27-41, 1998. https://doi.org/10.1109/83.650848
