# Changelog

## 0.2.2

### changed
- License changed to allow distribution on python package repositories

## 0.2.1

### added
- progress callback function can be supplied to `register_stack()` and `register_transform_stack()` functions via the `progress_callback` that is then called after every iteration (i.e. after each image registration)

### changed
- progress bar output not shown by default, has to be enabled by using the `verbose=True` parameter in the `register_stack()` and `register_transform_stack()` functions

## 0.2.0

### added
- bilinear transformation
