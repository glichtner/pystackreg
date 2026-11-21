Overview
==============

This section describes pystackreg's API.

The public interface consists of the class :class:`pystackreg.StackReg` and of
the functions in :mod:`pystackreg.util`. Everything else (in particular the
extension module ``pystackreg.turboreg`` and names starting with an
underscore) is an implementation detail and may change without notice.

.. currentmodule:: pystackreg

.. automodule:: pystackreg
   :no-undoc-members:

.. autosummary::

    StackReg
    util

Registration & Transformation
=============================

.. autoclass::  StackReg
   :members:

Utility functions
=================
.. currentmodule:: pystackreg.util
.. autofunction:: to_uint16
.. autofunction:: to_int_dtype
.. autofunction:: max_displacement
.. autofunction:: running_mean
.. autofunction:: simple_slice
