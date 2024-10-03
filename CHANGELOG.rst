Changelog
=========

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
