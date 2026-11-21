Changelog
=========

0.3.0 (unreleased)
------------------

This release modernizes the package. Apart from the corrections listed below,
results of registrations and transformations are unchanged: they are verified
bit for bit against 0.2.8.

License
.......
- pystackreg is now distributed under the **GNU General Public License,
  version 3 or later**, the license under which TurboReg and StackReg are
  distributed by their author. Versions up to 0.2.8 remain under the terms
  they were published with. See ``NOTICE``.

Added
.....
- Parallel processing: ``register_stack()``, ``transform_stack()`` and
  ``register_transform_stack()`` accept ``n_jobs`` (number of threads; ``-1``
  uses all CPUs). The native code releases the GIL, so threads are used and no
  data are copied between processes. Parallel results are identical to
  sequential results for all reference modes, including ``"previous"``
  (`issue #14 <https://github.com/glichtner/pystackreg/issues/14>`_,
  `PR #17 <https://github.com/glichtner/pystackreg/pull/17>`_,
  `PR #33 <https://github.com/glichtner/pystackreg/pull/33>`_).
- ``transform_stack()`` and ``register_transform_stack()`` accept ``out``: a
  preallocated array of another dtype (integer dtypes are rounded and clipped),
  a ``numpy.memmap``, or the input itself for in-place operation
  (`issue #28 <https://github.com/glichtner/pystackreg/issues/28>`_).
- ``StackReg.transform_points()`` applies a transformation to point
  coordinates, also for the bilinear transformation
  (`issue #42 <https://github.com/glichtner/pystackreg/issues/42>`_).
- ``pystackreg.util.max_displacement()`` measures the amplitude of
  transformations, e.g. to detect implausible registration results
  (`issue #41 <https://github.com/glichtner/pystackreg/issues/41>`_).
- ``transform_stack()`` accepts ``verbose`` and ``progress_callback``;
  ``register_transform_stack()`` accepts ``suppress_axis_warning``.
- Type annotations (``py.typed``), support for free-threaded CPython builds.
- Validation suite; the documentation reports a comparison with the original
  TurboReg.

Changed
.......
- Python 3.10 or newer is required; Python 3.10 - 3.15 are tested. Wheels are
  built against NumPy 2 and run with NumPy >= 1.22.
- Memory use of stack functions is much lower
  (`issue #28 <https://github.com/glichtner/pystackreg/issues/28>`_):
  the time-axis check, the moving average and ``transform_stack()`` no longer
  create temporary copies of the whole stack.
- ``transform()`` is about three times faster (1.6 times for the bilinear
  transformation); it no longer computes image pyramids and gradients that
  only the registration needs. The speed of ``register()`` is unchanged.
- The check that warns about a wrong ``axis`` is evaluated in chunks.
- Binary wheels are no longer built for 32-bit platforms (Linux i686,
  Windows win32); wheels for Linux aarch64 and macOS x86-64 are added.
- Progress bars use ``tqdm.auto`` (terminal or notebook widget,
  `PR #11 <https://github.com/glichtner/pystackreg/pull/11>`_), are labelled,
  and are closed when an error occurs. ``tqdm`` is imported only when
  ``verbose=True``. ``register_transform_stack(verbose=True)`` shows a second
  bar for the transformation.
- Invalid arguments raise ``ValueError``, ``TypeError`` or ``RuntimeError``
  with specific messages instead of a bare ``Exception`` (all three are
  subclasses of ``Exception``).
- Images containing NaN or infinite values are rejected by ``register()`` and
  ``register_stack()`` with a ``ValueError``. Previously such images silently
  produced an identity transformation.
- Images smaller than 3 x 3 pixels, empty stacks, non-finite transformation
  matrices and matrices of the wrong shape are rejected with a ``ValueError``.
  Previously these cases read out-of-bounds memory or crashed.
- The moving average (``moving_average`` > 1) of floating-point stacks is
  computed in double precision, so results for such stacks can differ in the
  last digits from 0.2.8. Integer stacks are unaffected.
- Compiled with ``-ffp-contract=off``: results no longer depend on whether
  the CPU architecture fuses multiply-add operations (affects arm64 builds in
  the last bits).
- Project metadata moved to ``pyproject.toml``; development uses Ruff, mypy
  and pytest.

Fixed
.....
- ``reference="previous"`` combined with ``moving_average`` > 1 registered
  the first frame against the *last* frame of the stack, which shifted the
  whole registered stack.
- ``moving_average`` > 1 gave wrong results for ``axis`` other than 0
  (``util.running_mean()`` ignored the axis when taking differences).
- Memory leak in the extension module: every call to ``register()`` leaked
  the three returned arrays (about 430 bytes).
- Crash (NULL pointer dereference) when an image could not be converted to
  floating point, e.g. for complex data.
- Out-of-bounds reads for images with a side of one or two pixels.
- Undefined behavior (float-to-int conversion of NaN/out-of-range values) in
  the registration of constant images and in ``transform()`` with extreme
  matrices; conversions now follow the semantics of the Java original.
- ``get_matrix()`` before a registration raised an unrelated ``TypeError``.
- Images with a side shorter than 24 pixels were not registered reliably with
  the scaled-rotation and affine transformations: the registration used a
  value that was never set, so the result depended on the contents of memory
  and was usually the identity. Such images are now registered. This differs
  from TurboReg, which returns the identity in this case; results for larger
  images are unchanged.

Removed
.......
- Support for Python < 3.10 and the Python 2 code paths of the extension.
- ``pystackreg/conftest.py`` is no longer installed with the package.

0.2.8
-----

Fixed
.....
- Add NumPy 2 support

0.2.7
-----

Fixed
.....
- Axis argument not used for method "mean" in register_stack() (`PR #26 <https://github.com/glichtner/pystackreg/pull/26>`_)

0.2.6
-----

Added
.....
- Exposing `simple_slice` and `running_mean` functions in the util package
- Added conversion function to any integer dtype

0.2.5
-----

Fixed
.....
- Compilation in environments without NumPy

0.2.3
-----

Added
.....
- Added example data and tutorial notebook
- Added unit tests
- Additional documentation
- Detection of time series axis in stacks – will raise a warning if supplied axis in stack registration does not correspond to the detected axis

Changed
.....~~
- `progress_callback` function now gets called with the iteration number, not the iteration index (iteration number = iteration index + 1)

Fixed
.....
- Fixed exception when using a different axis than 0 for registering stacks

0.2.2
-----

Changed
.....~~
- License changed to allow distribution on Python package repositories

0.2.1
-----

Added
.....
- Progress callback function can be supplied to `register_stack()` and `register_transform_stack()` functions via the `progress_callback` parameter. It is called after every iteration (i.e., after each image registration).

Changed
.....~~
- Progress bar output is not shown by default, has to be enabled by using the `verbose=True` parameter in the `register_stack()` and `register_transform_stack()` functions

0.2.0
-----

Added
.....
- Bilinear transformation
