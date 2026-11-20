Installation and troubleshooting
================================

Supported platforms
-------------------

.. list-table::
   :header-rows: 1

   * - Component
     - Supported
   * - Python
     - 3.10 - 3.15 (CPython), including the free-threaded build of 3.14
   * - NumPy
     - 1.22 or newer, including NumPy 2
   * - Operating systems (binary wheels)
     - Linux x86-64 and aarch64 (glibc and musl), macOS x86-64 and arm64,
       Windows x86-64
   * - Other platforms
     - installation from source with a C++11 compiler

Installing
----------

.. code-block:: console

   pip install pystackreg

or, with conda/mamba:

.. code-block:: console

   conda install -c conda-forge pystackreg

From source (requires a C++ compiler; NumPy is installed automatically as a
build dependency):

.. code-block:: console

   pip install git+https://github.com/glichtner/pystackreg

Check the installation with

.. code-block:: console

   python -c "import pystackreg; print(pystackreg.__version__)"

Installation problems
---------------------

``error: command 'gcc' failed`` / ``Microsoft Visual C++ 14.0 or greater is required``
    pip found no binary wheel for your platform or Python version and tried to
    compile the package. Upgrade pip (``python -m pip install -U pip``) so
    that current wheel tags are recognized, or install a compiler (Linux:
    ``build-essential`` and the Python headers, e.g. ``python3-dev``; macOS:
    ``xcode-select --install``; Windows: "Build Tools for Visual Studio").

``ImportError: cannot import name 'turboreg' from partially initialized module 'pystackreg'``
    Python imported the *source directory* named ``pystackreg`` instead of the
    installed package, typically because the interpreter was started in a
    checkout of the repository. Change to another directory, or build the
    extension in place with ``pip install -e .``.

``ImportError: numpy.core.multiarray failed to import`` or a message about a module compiled with NumPy 1.x
    An old pystackreg build (0.2.7 or earlier) is used with NumPy 2. Upgrade
    pystackreg.

Errors and warnings at run time
-------------------------------

``UserWarning: Detected axis … as the possible time axis``
    The stack varies least along another axis than the one given by ``axis``.
    Usually the frames are not on the first axis: pass ``axis=…``. For short
    stacks, or stacks that change a lot over time, the heuristic can be
    wrong; pass ``suppress_axis_warning=True`` then.

``ValueError: … contains NaN or infinite values; such images cannot be registered``
    Replace or mask these pixels first, e.g. ``np.nan_to_num(img, nan=np.nanmedian(img))``.

``ValueError: Bilinear stack transformation not supported with reference == "previous"``
    Bilinear transformations cannot be chained. Use ``reference="first"`` or
    ``"mean"``.

``ValueError: Number of saved transformation matrices does not match stack length``
    ``transform_stack`` needs one matrix per frame of the stack that is being
    transformed, along the given ``axis``.

``RuntimeError: No transformation available`` / ``No transformation matrices given``
    Call ``register``/``register_stack`` first, or pass the matrices.

The registered images contain negative values
    Expected: cubic interpolation overshoots at sharp edges. See
    :ref:`limitations`.

The result is not aligned, or worse than before
    Check, in this order: that ``axis`` is the frame axis; that the
    displacement is within the capture range of the method (see
    :ref:`limitations`); that both images show the same contrast; and that
    the chosen motion model is not more flexible than necessary. For
    ``reference="first"`` on long recordings, try ``"previous"`` or
    ``"mean"``. :func:`pystackreg.util.max_displacement` helps to find the
    frames where the registration went wrong.

The transformed image is shifted in the wrong direction when I apply the matrix to points
    The matrix maps reference coordinates onto coordinates of the moving
    image, and coordinates are (x, y), not (row, column). Use
    :meth:`pystackreg.StackReg.transform_points`; see :ref:`matrix-convention`.

Not enough memory
    See :doc:`performance`: pass ``out=`` to ``transform_stack`` and reduce
    ``n_jobs``.

Results differ slightly from ImageJ
    ImageJ shows and saves results of StackReg in the bit depth of the input.
    Compare after rounding (ImageJ rounds half up), or compare the 32-bit
    output of TurboReg. See :doc:`validation` for what agreement to expect.
    Also make sure that TurboReg was run in "accurate" mode.

Getting help
------------

Please open an issue at https://github.com/glichtner/pystackreg/issues and
include the version numbers of pystackreg, NumPy and Python, the operating
system, and a minimal example.
