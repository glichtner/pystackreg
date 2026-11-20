Validation and reproducibility
==============================

pystackreg is a port of existing software. Its value depends on reproducing
the original faithfully, so this is tested at three levels. All tests are part
of the test suite (``tests/``) and run in continuous integration.

.. contents::
   :local:
   :depth: 1

Recovery of known transformations
---------------------------------

``tests/test_accuracy.py`` uses analytic test scenes (sums of Gaussian blobs,
``tests/helpers.py``). A transformed copy of such a scene is obtained by
evaluating the function at transformed coordinates, so the true transformation
is known exactly and no resampling error enters the test data. Accuracy is
measured as the largest distance between the positions to which the estimated
and the true transformation map the four image corners.

.. list-table:: Corner error in pixels (pystackreg 0.3.0, Linux x86-64)
   :header-rows: 1

   * - Image size (rows x columns)
     - translation
     - rigid body
     - scaled rotation
     - affine
     - bilinear
   * - 96 x 128
     - < 0.0001
     - < 0.0001
     - 0.0005
     - < 0.0001
     - 0.0001
   * - 128 x 96
     - < 0.0001
     - < 0.0001
     - 0.0004
     - < 0.0001
     - 0.0001
   * - 150 x 110
     - < 0.0001
     - < 0.0001
     - 0.0002
     - < 0.0001
     - 0.0001
   * - 64 x 64
     - 0.0001
     - 0.0001
     - 0.0002
     - 0.0003
     - 0.0010

Further properties covered by the suite: identical images give the identity;
sub-pixel shifts between 0.05 and 0.75 pixel are recovered to better than
0.001 pixel; swapping the two images gives the inverse transformation; the
result does not depend on the intensity scale or offset; the matrices have the
structure of their model (e.g. a pure rotation for the rigid-body model); the
chained transformations of a drifting stack match the ground truth for all
reference modes.

These numbers describe noise-free, well-textured synthetic images. Accuracy
on real data depends on noise, image content and on how well the motion model
fits.

Comparison with the original TurboReg
-------------------------------------

For release 0.3.0 the unmodified TurboReg plugin (version "November 30,
2024") was run headlessly with ImageJ 1.54p on 20 image pairs: analytic scenes
with 32-bit and 16-bit pixels and two frames of a real time-lapse recording,
for all five transformation types. TurboReg was called in the two steps in
which StackReg calls it: alignment of the images cropped by one row and
column, then transformation of the full image. Its refined landmarks and
transformed images were compared with those of pystackreg.

.. list-table:: Largest difference to TurboReg over all image pairs
   :header-rows: 1

   * - Model
     - Refined landmarks
     - Transformed image
   * - scaled rotation, affine, bilinear
     - bit-identical
     - bit-identical
   * - translation
     - 4 x 10\ :sup:`-11` pixel
     - bit-identical
   * - rigid body
     - 6 x 10\ :sup:`-8` pixel
     - 6 x 10\ :sup:`-8` relative (last bit of single precision)

Two known causes explain the differences for the translation and rigid-body
models:

* TurboReg stores images and gradients in single precision. In the Hessian of
  these two models it multiplies two such values in single precision
  (``float * float`` in Java, twelve statements). The C++ port holds the same
  single-precision values in double-precision variables and forms the products
  in double precision.
* The rigid-body model uses ``sin``, ``cos`` and ``atan2``, whose last bit
  differs between the Java runtime and the C library.

Both effects are far below any practical relevance, and they are left as they
are so that pystackreg keeps reproducing the results of its earlier releases.

This comparison needs a Java runtime and the TurboReg sources and is not part
of the test suite. The regression tests described below tie later versions to
the results for which it was made.

The test suite itself contains a coarser comparison:
``tests/test_imagej_reference.py`` (part of pystackreg since 0.2.3) compares
complete registered stacks with stacks that were registered with the StackReg
and TurboReg plugins in ImageJ and saved as 16-bit images; they agree within
one grey level, the rounding uncertainty of the stored data.

Independent check of the interpolation
--------------------------------------

``tests/test_transform.py`` compares :meth:`~pystackreg.StackReg.transform`
with SciPy's cubic B-spline interpolation (``scipy.ndimage.map_coordinates``,
``order=3``, ``mode="reflect"``), an implementation that shares no code with
TurboReg. Wherever TurboReg produces data the two agree to single precision.

Stability across pystackreg versions
------------------------------------

``tests/data/baseline_v0.2.8.npz`` holds 76 result arrays (matrices, landmarks
and images) computed with the unmodified release 0.2.8. On the platform where
the baseline was generated (Linux x86-64, GCC) the current code reproduces
all of them **bit for bit**, although the native code was refactored in 0.3.0.
Set ``PYSTACKREG_EXACT_BASELINE=1`` to make the tests insist on this;
otherwise they allow 10\ :sup:`-6` pixel for landmarks and matrices and
10\ :sup:`-5` relative for image intensities, to accommodate other compilers
and math libraries.

Two deliberate exceptions are documented in the changelog: the registration
of floating-point stacks with ``moving_average`` > 1 (now averaged in double
precision), and the cases that were defects in 0.2.8.

If results ever have to change deliberately, the reason is documented in the
changelog and a new baseline is generated with ``tests/generate_baseline.py``
from a trusted build; the previous baseline is kept for reference.

Reproducibility
---------------

* The algorithm is deterministic: there is no random initialization, and the
  number of iterations is bounded.
* Parallel execution (``n_jobs`` > 1) gives bit-identical results to
  sequential execution, for every reference mode and independent of the
  number of threads and of scheduling. Frames are independent computations;
  for ``reference="previous"`` the pairwise results are chained afterwards
  in frame order.
* The extension is compiled without value-changing optimizations (no
  ``-ffast-math``) and with ``-ffp-contract=off``, so that results do not
  depend on whether the CPU has fused multiply-add instructions. Remaining
  platform dependence comes from the C math library (``sin``, ``cos``,
  ``atan2``, ``pow``) and is of the order of the last bits. It can be
  amplified to roughly 10\ :sup:`-7` pixel by the iterative optimization;
  bit-identical results across operating systems are therefore not
  guaranteed.
* To reproduce an analysis exactly, record the version of pystackreg and
  store the transformation matrices together with the results; applying
  stored matrices does not involve the optimizer.
* Input data are converted to double precision. Integer data of any width up
  to 32 bit and single-precision data are represented exactly.

Memory safety
-------------

The native code is exercised under AddressSanitizer and
UndefinedBehaviorSanitizer (including float-to-int conversion checks) with a
battery of edge cases - minimal image sizes, constant images, non-finite
values, degenerate matrices, non-contiguous and read-only arrays. Version
0.3.0 is clean under these tools; the defects they revealed in 0.2.8 are
listed in the changelog. ``tools/sanitize.sh`` builds the extension with both
sanitizers and runs the test suite against it.
