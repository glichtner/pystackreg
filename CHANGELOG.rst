#########
Changelog
#########

0.2.5
=====

fixed
-----
- Compilation in environment without numpy

0.2.3
=====

added
-----
- Added example data and tutorial notebook
- Added unit tests
- Additional documentation
- Detection of time series axis in stacks - will raise a warning if supplied axis in
  stack registration does not correspond to the detected axis

changed
-------
- `progress_callback` function now gets called with the iteration number, not
  iteration index (iteration number = iteration index + 1)

fixed
-----
- Fixed exception when using a different axis than 0 for registering stacks

0.2.2
=====

changed
-------
- License changed to allow distribution on python package repositories

0.2.1
=====

added
-----
- progress callback function can be supplied to `register_stack()` and `register_transform_stack()` functions via the `progress_callback` that is then called after every iteration (i.e. after each image registration)

changed
-------
- progress bar output not shown by default, has to be enabled by using the `verbose=True` parameter in the `register_stack()` and `register_transform_stack()` functions

0.2.0
=====

added
-----
- bilinear transformation
