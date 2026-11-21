Performance and large stacks
============================

This page explains how to process stacks faster and with less memory, and what
to expect. The numbers were measured on a laptop (Intel Core i7-1260P) and are
meant for orientation; measure with your own data and hardware.

Parallel processing
-------------------

:meth:`~pystackreg.StackReg.register_stack`,
:meth:`~pystackreg.StackReg.transform_stack` and
:meth:`~pystackreg.StackReg.register_transform_stack` take an ``n_jobs``
argument::

    sr = StackReg(StackReg.RIGID_BODY)
    registered = sr.register_transform_stack(stack, reference="previous", n_jobs=-1)

``n_jobs=1`` (the default) processes the frames one after the other in the
calling thread, exactly as earlier versions did. A larger value uses that many
threads, ``-1`` uses all available CPU cores, ``-2`` all but one.

How it works
............

The registration and resampling code is native C++ that keeps all of its state
in objects local to one call. pystackreg releases Python's global interpreter
lock while it runs, so several frames are processed truly in parallel by
ordinary threads. Threads share the memory of the stack: nothing is copied or
serialized, in contrast to process-based parallelism.

Every frame is an independent task:

* ``reference="first"`` and ``"mean"``: each frame is registered to the same,
  fixed reference image.
* ``reference="previous"``: frame *i* is registered to frame *i - 1* of the
  *unregistered* stack. These pairwise registrations do not depend on each
  other. After all of them have finished, the pairwise matrices are multiplied
  in frame order, :math:`T_i = P_i\,T_{i-1}`, exactly as in the sequential
  computation.
* :meth:`~pystackreg.StackReg.transform_stack`: each frame is resampled with
  its own matrix.

The result therefore does **not** depend on ``n_jobs``, on the number of CPU
cores or on the order in which the threads finish: it is bit-identical to the
sequential result. The test suite checks this for all transformation types and
reference modes.

If a frame raises an error (for example because it contains NaN), frames that
have not been started are cancelled, the error is raised in the calling
thread, and the :class:`~pystackreg.StackReg` object remains usable. The same
happens on ``KeyboardInterrupt``; a frame that is already being processed is
finished first, which takes at most the time of one frame.

What to expect
..............

The speed-up depends on the machine. On the laptop CPU used for these
measurements (4 performance cores and 8 efficiency cores, 16 hardware threads),
registering 64 frames of 512 x 512 pixels took 4.2 s with one thread. With 16
threads it took 0.76 s (5.5 times faster) as long as the processor was not
throttled, and about 1.8 s (2.3 times faster) under sustained load, when the
power management of the laptop slows all cores down; transformation scaled
similarly. Apart from such throttling, scaling is less than linear because the
cores are not equally fast, because the clock frequency drops when all cores
are busy, and because the algorithm is limited by memory bandwidth. These
numbers characterize that laptop more than they characterize pystackreg: no
workstation or server processor was available for the measurements, so measure
with your own data and hardware.

Compared with a pool of worker *processes*, threads registered the same stack
in about the same time (1.20 s against 1.33 s with 8 workers) but needed no
worker start-up and copied no pixel data, whereas the process pool had to
serialize 158 MB.

Each thread needs working memory for one registration, about 14 images of the
frame size in double precision (28 MB for 512 x 512 pixels). For very large
frames, choose ``n_jobs`` accordingly.

Using pystackreg from your own threads
......................................

Because the lock is released, pystackreg can also be called from threads that
you manage yourself (``concurrent.futures.ThreadPoolExecutor``, Dask's threaded
scheduler). Use a separate :class:`~pystackreg.StackReg` object per thread:
an object stores the result of its last registration.

Memory
------

A registration needs little memory beyond the stack itself. The main consumer
is the *result* of the transformation: by default
:meth:`~pystackreg.StackReg.transform_stack` returns a ``float64`` array, which
is four times as large as a 16-bit stack and eight times as large as an 8-bit
stack.

.. list-table:: Peak memory in addition to the input (200 frames of
                512 x 512 pixels, 16 bit, 100 MB)
   :header-rows: 1

   * - Operation
     - 0.2.8
     - 0.3.0
   * - ``register_stack``
     - 400 MB
     - 34 MB
   * - ``register_stack(moving_average=5)``
     - 1210 MB
     - 34 MB
   * - ``register_stack(n_jobs=8)``
     - \-
     - 228 MB
   * - ``transform_stack`` (float64 result)
     - 500 MB
     - 406 MB
   * - ``transform_stack(out=...)`` (float32, uint16 or in place)
     - \-
     - 6 MB

The ``out`` argument
....................

``out`` makes :meth:`~pystackreg.StackReg.transform_stack` (and
:meth:`~pystackreg.StackReg.register_transform_stack`) write into an existing
array instead of allocating the ``float64`` result:

.. code-block:: python

    # single precision: half the size, no loss (the result of TurboReg has
    # single precision anyway)
    out = np.empty(stack.shape, dtype=np.float32)
    sr.transform_stack(stack, out=out)

    # same integer type as the input: values are rounded and clipped
    out = np.empty_like(stack)
    sr.transform_stack(stack, out=out)

    # in place: no second stack at all
    sr.transform_stack(stack, out=stack)

For integer arrays the values are rounded to the nearest integer and clipped
to the range of the type, as :func:`pystackreg.util.to_int_dtype` does.
Negative overshoot of the interpolation is thus set to 0 for unsigned types.

Stacks larger than memory
.........................

Input and output can be memory-mapped files; pystackreg then reads and writes
one frame at a time (per thread):

.. code-block:: python

    import numpy as np
    from pystackreg import StackReg

    stack = np.load("recording.npy", mmap_mode="r")          # on disk
    out = np.lib.format.open_memmap(
        "registered.npy", mode="w+", dtype=np.float32, shape=stack.shape
    )

    sr = StackReg(StackReg.RIGID_BODY)
    tmats = sr.register_stack(stack, reference="previous", n_jobs=4, verbose=True)
    np.save("tmats.npy", tmats)
    sr.transform_stack(stack, out=out, n_jobs=4, verbose=True)
    out.flush()

TIFF files can be memory-mapped with ``tifffile.memmap`` if they are
uncompressed and contiguous, Zarr or HDF5 data can be processed in blocks:
register a block, keep the last frame as the reference for the next block
(``reference="previous"``), or register all blocks to one reference image
(``reference="first"`` with the reference prepended).

Notes:

* ``reference="mean"`` reads the whole stack once to compute the mean image.
* The time-axis check reads the whole stack once as well. For stacks on slow
  storage, pass ``suppress_axis_warning=True`` to skip it.
* ``moving_average=N`` reads ``N`` frames per registered frame.

Progress reporting
------------------

``verbose=True`` shows a `tqdm <https://tqdm.github.io>`_ progress bar: a
text bar in a terminal and a widget in Jupyter (``tqdm.auto``).
``register_transform_stack`` shows one bar for the registration and one for
the transformation. Nothing is printed by default.

For graphical user interfaces or logging, pass a function as
``progress_callback``. It is called after every frame with the keyword
arguments ``current_iteration`` (number of frames finished so far) and
``end_iteration`` (total number of frames to process)::

    def report(current_iteration, end_iteration):
        print(f"{current_iteration} of {end_iteration} frames")

    sr.register_stack(stack, progress_callback=report, n_jobs=4)

The callback always runs in the thread that called pystackreg, also when
several worker threads are used, so it may update a user interface directly.
In parallel runs the count increases in the order in which frames finish. An
exception raised by the callback stops the processing.

Speed of a single registration
------------------------------

Run time is proportional to the number of pixels. On one core of that
laptop, registering a pair of 512 x 512 images takes about 67 ms and
resampling an image about 12 ms; the transformation type matters little
(translation is fastest).

Version 0.3.0 did not change the speed of a registration; resampling
(``transform``) became about three times faster (1.6 times for the bilinear
transformation) because it no longer computes image pyramids and gradients
that only the registration needs.

Why there is no GPU backend
---------------------------

The optimizer is an iterative procedure with data-dependent control flow that
runs on an image pyramid; a single frame takes tens of milliseconds on one CPU
core. A GPU implementation would be a separate implementation of the algorithm
with different rounding behaviour, which would give up the agreement with
TurboReg that pystackreg is validated for. Throughput for stacks is instead
increased with CPU threads, which requires no additional dependency.
