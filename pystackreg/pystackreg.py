# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2018-2026 Gregor Lichtner and PyStackReg contributors
"""Image and image-stack registration: the :class:`StackReg` class."""

from __future__ import annotations

import operator
import warnings
from collections.abc import Callable
from typing import Any

import numpy as np
import numpy.typing as npt

from pystackreg import turboreg  # type: ignore[attr-defined]
from pystackreg._internal import (
    FrameSource,
    Progress,
    check_output_array,
    detect_time_axis,
    resolve_n_jobs,
    run_tasks,
    store_frame,
)
from pystackreg.util import simple_slice

__all__ = ["StackReg"]

#: Smallest image side (in pixels) that can be registered.
_MIN_REGISTRATION_SIZE = 3


class StackReg:
    """
    Python implementation of the ImageJ/Fiji StackReg plugin
    (https://bigwww.epfl.ch/thevenaz/stackreg/ )

    One instance holds one transformation type and the result of the most
    recent registration. Instances are cheap to create; they are not meant to
    be shared between threads that register concurrently (use one instance per
    thread, or the ``n_jobs`` argument of the stack functions).

    :type transformation: int
    :param transformation: One of :attr:`TRANSLATION`, :attr:`RIGID_BODY`,
                           :attr:`SCALED_ROTATION`, :attr:`AFFINE`,
                           :attr:`BILINEAR`
    """

    TRANSLATION = 2
    """
    Translation
    """

    RIGID_BODY = 3
    """
    Rigid body (translation + rotation)
    """

    SCALED_ROTATION = 4
    """
    Scaled rotation (translation + rotation + scaling)
    """

    AFFINE = 6
    """
    Affine (translation + rotation + scaling + shearing)
    """

    BILINEAR = 8
    """
    Bilinear (non-linear transformation; does not preserve straight lines)
    """

    _valid_transformations = (
        TRANSLATION,
        RIGID_BODY,
        SCALED_ROTATION,
        AFFINE,
        BILINEAR,
    )

    _valid_references = ("previous", "first", "mean")

    _is_registered = False

    def __init__(self, transformation: int) -> None:
        if transformation not in self._valid_transformations:
            raise ValueError(
                f"Invalid transformation {transformation!r}; use one of "
                "StackReg.TRANSLATION, RIGID_BODY, SCALED_ROTATION, AFFINE, "
                "BILINEAR"
            )

        self._transformation = transformation
        self._m: np.ndarray | None = None
        self._tmats: np.ndarray | None = None
        self._refpts: np.ndarray | None = None
        self._movpts: np.ndarray | None = None

    def is_registered(self) -> bool:
        """
        Indicates whether register() was already called and a transformation matrix
        was calculated

        :rtype:  bool
        :return: True if a transformation matrix was already calculated
        """
        return self._is_registered

    def register(self, ref: npt.ArrayLike, mov: npt.ArrayLike) -> np.ndarray:
        """
        Registers an image to a reference image: Only the transformation matrix
        will be calculated, the image will not be transformed
        (use transform() or register_transform() ).

        As in StackReg, the last row and the last column of both images are
        not used for the registration.

        :type ref: array_like (M, N)
        :param ref: Reference image (static)

        :type mov: array_like (M, N)
        :param mov:
            The image that should be aligned to the reference image. Must have
            the same shape as ``ref``.

        :rtype:  ndarray(3,3) or ndarray(4,4) for bilinear transformation
        :return: Transformation matrix. It maps coordinates ``(x, y)``
                 (column, row) of the reference image onto coordinates of
                 ``mov``; see :meth:`transform_points`.

        :raises ValueError: if the images are not 2-D, differ in shape, are
            smaller than 3 x 3 pixels or contain NaN or infinite values
        """
        self._m, self._refpts, self._movpts = self._register_pair(
            np.asarray(ref), np.asarray(mov)
        )
        self._is_registered = True

        return self.get_matrix()

    def _register_pair(
        self, ref: np.ndarray, mov: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Validate and register one pair of images without touching the state of
        this object (which makes it safe to call from several threads).

        :return: (compact TurboReg matrix, reference points, moving points)
        """
        for name, image in (("ref", ref), ("mov", mov)):
            if image.ndim != 2:
                raise ValueError(
                    f"{name} must be a two-dimensional image, got an array "
                    f"with {image.ndim} dimension(s)"
                )
        if ref.shape != mov.shape:
            raise ValueError(
                f"ref and mov must have the same shape, got {ref.shape} and {mov.shape}"
            )
        if min(ref.shape) < _MIN_REGISTRATION_SIZE:
            raise ValueError(
                f"Images must be at least {_MIN_REGISTRATION_SIZE} x "
                f"{_MIN_REGISTRATION_SIZE} pixels for registration, got "
                f"{ref.shape[0]} x {ref.shape[1]}"
            )
        for name, image in (("ref", ref), ("mov", mov)):
            if image.dtype.kind == "f" and not np.isfinite(image).all():
                raise ValueError(
                    f"{name} contains NaN or infinite values; such images "
                    "cannot be registered"
                )

        # StackReg hands TurboReg the images without their last row and column
        return turboreg._register(ref[:-1, :-1], mov[:-1, :-1], self._transformation)

    def transform(
        self, mov: npt.ArrayLike, tmat: npt.ArrayLike | None = None
    ) -> np.ndarray:
        """
        Transform an image according to a previous registration.

        Either a transformation matrix has to be explicitly supplied or
        register() has to be called before calling transform().

        The image is resampled with cubic B-spline interpolation. Positions
        that map outside of ``mov`` are set to 0. Because of the cubic
        interpolation, values can slightly over- or undershoot the range of
        the input (e.g. become negative next to sharp edges); see
        :func:`pystackreg.util.to_uint16` for a conversion that clips.

        :type mov: array_like (M, N)
        :param mov: The image that will be transformed

        :type tmat:  ndarray(3,3) or ndarray(4,4) for bilinear transformation,
                     optional
        :param tmat: The transformation matrix

        :rtype:  ndarray (M, N), float64
        :return: Transformed image - will be of the same shape as the input image
                 (cropping may occur)
        """
        if tmat is None:
            if self._m is None:
                raise RuntimeError(
                    "No transformation available: call register() first or "
                    "pass a transformation matrix"
                )
            m = self._m
        else:
            m = self._matrix_long_to_short(np.asarray(tmat))

        return turboreg._transform(mov, m)

    def register_transform(self, ref: npt.ArrayLike, mov: npt.ArrayLike) -> np.ndarray:
        """
        Register and transform an image to a reference image.

        :type ref: array_like (M, N)
        :param ref: Reference image (static)

        :type mov: array_like (M, N)
        :param mov: The image that should be aligned to the reference image

        :rtype:  ndarray (M, N), float64
        :return: Transformed image - will be of the same shape as the input image
                 (cropping may occur)
        """
        self.register(ref, mov)
        return self.transform(mov)

    def get_matrix(self) -> np.ndarray:
        """
        Get the current transformation matrix

        :rtype:  ndarray(3,3) or ndarray(4,4) for bilinear transformation
        :return: The transformation matrix
        """
        if self._m is None:
            raise RuntimeError(
                "No transformation available: call register() or set_matrix() first"
            )
        return self._matrix_short_to_long(self._m)

    def set_matrix(self, mat: npt.ArrayLike) -> None:
        """
        Sets the current transformation matrix

        :type mat: ndarray(3,3) or ndarray(4,4) for bilinear transformation
        :param mat: The transformation matrix
        """
        mat = np.asarray(mat)
        exp_shape = (4, 4) if self._transformation == self.BILINEAR else (3, 3)

        if exp_shape != mat.shape:
            raise ValueError(
                f"Invalid shape of transformation matrix: Expected {exp_shape}, "
                f"got {mat.shape}"
            )

        self._m = self._matrix_long_to_short(mat)

        self._is_registered = True

    def _matrix_short_to_long(self, m: np.ndarray) -> np.ndarray:
        """
        Converts the transformation matrix from the short form used by TurboReg to
        the canonical form from linear algebra.

        :type m: array_like
        :param m: TurboReg transformation matrix

        :rtype:  ndarray(3,3)
        :return: Canonical transformation matrix
        """

        if self._transformation == self.TRANSLATION:
            mat = np.identity(3).astype(np.double)
            mat[0:2, 2] = m[:, 0]
        elif self._transformation in [
            self.RIGID_BODY,
            self.SCALED_ROTATION,
            self.AFFINE,
        ]:
            mat = np.identity(3).astype(np.double)
            mat[0:2, :] = m[:, [1, 2, 0]]
        else:  # BILINEAR
            mat = np.identity(4).astype(np.double)
            mat[0:2, :] = m[:, [1, 2, 3, 0]]

        return mat

    def _matrix_long_to_short(self, mat: np.ndarray) -> np.ndarray:
        """
        Converts the transformation matrix from the canonical form from linear algebra
        to the short form used by TurboReg.

        :type mat: ndarray(3,3)
        :param mat: Canonical transformation matrix

        :rtype:  array_like
        :return: TurboReg transformation matrix
        """
        n_cols = 4 if self._transformation == self.BILINEAR else 3
        if mat.ndim != 2 or mat.shape[0] < 2 or mat.shape[1] != n_cols:
            raise ValueError(
                f"Invalid shape of transformation matrix: Expected "
                f"({n_cols}, {n_cols}), got {mat.shape}"
            )

        if self._transformation == self.TRANSLATION:
            m = mat[0:2, 2].reshape((2, 1))
        elif self._transformation in [
            self.RIGID_BODY,
            self.SCALED_ROTATION,
            self.AFFINE,
        ]:
            m = mat[0:2, [2, 0, 1]]
        else:  # BILINEAR
            m = mat[0:2, [3, 0, 1, 2]]

        return m.astype(np.double)

    @staticmethod
    def _detect_time_axis(img: np.ndarray) -> int:
        """
        Detects the time axis of a movie/stack by returning the dimension with the
        lowest average variability

        :param img: 3-D stack
        :return: index of the axis with the lowest average variance
        """
        return detect_time_axis(img)

    def get_points(self) -> tuple[np.ndarray | None, np.ndarray | None]:
        """
        Returns the pairs of corresponding points from which the
        transformation matrix can be calculated.

        :rtype:  (ndarray, ndarray)
        :return: (Point coordinates of reference image, Point coordinates of image to
                 be aligned)
        """
        return self._refpts, self._movpts

    def transform_points(
        self,
        points: npt.ArrayLike,
        tmat: npt.ArrayLike | None = None,
        *,
        inverse: bool = False,
    ) -> np.ndarray:
        """
        Apply a transformation to point coordinates instead of to an image.

        By default the points are moved in the same way as :meth:`transform`
        moves image content: a feature located at ``points`` in the
        unregistered image ``mov`` is found at the returned coordinates in the
        registered image (i.e. in the frame of the reference image).

        Coordinates are given as ``(x, y)`` = ``(column, row)``, with the
        origin in the centre of the top-left pixel. Use ``points[..., ::-1]``
        to convert from/to NumPy's ``(row, column)`` index order.

        .. note:: The transformation matrix itself maps *reference*
           coordinates onto *moving-image* coordinates (this is the direction
           needed to resample an image). Multiplying a point by the matrix
           therefore corresponds to ``inverse=True``.

        :type points: array_like (..., 2)
        :param points: Point coordinates ``(x, y)``

        :type tmat:  ndarray(3,3) or ndarray(4,4) for bilinear transformation,
                     optional
        :param tmat: The transformation matrix; defaults to the matrix of the
                     last registration

        :type inverse: bool, optional
        :param inverse:
            If False (default), map coordinates of the moving image to the
            registered/reference frame. If True, map coordinates of the
            reference frame to the moving image.

        :rtype:  ndarray (..., 2), float64
        :return: Transformed point coordinates ``(x, y)``
        """
        pts = np.asarray(points, dtype=np.float64)
        if pts.ndim == 0 or pts.shape[-1] != 2:
            raise ValueError(
                "points must have shape (..., 2) holding (x, y) coordinates, "
                f"got {pts.shape}"
            )
        if tmat is None:
            mat = self.get_matrix()
        else:
            # round trip through the compact form: use exactly the part of
            # the matrix that transform() uses
            mat = self._matrix_short_to_long(
                self._matrix_long_to_short(np.asarray(tmat))
            )

        x, y = pts[..., 0], pts[..., 1]
        if self._transformation == self.BILINEAR:
            u, v = _bilinear_map(mat, x, y) if inverse else _bilinear_solve(mat, x, y)
        else:
            if not inverse:
                mat = np.linalg.inv(mat)
            u = mat[0, 0] * x + mat[0, 1] * y + mat[0, 2]
            v = mat[1, 0] * x + mat[1, 1] * y + mat[1, 2]
        return np.stack([u, v], axis=-1)

    def register_stack(
        self,
        img: npt.ArrayLike,
        reference: str = "previous",
        n_frames: int = 1,
        axis: int = 0,
        moving_average: int = 1,
        verbose: bool = False,
        progress_callback: Callable[..., Any] | None = None,
        suppress_axis_warning: bool = False,
        *,
        n_jobs: int | None = 1,
    ) -> np.ndarray:
        """
        Register a stack of images (movie).
        Note that this function will not transform the image but only calculate the
        transformation matrices.
        For transformation use transform_stack() after this function or use
        register_transform_stack() for a single call.

        :type img: array_like(Ni..., Nj..., Nk...)
        :param img: The 3D stack of images that should be aligned

        :type reference: string, one of ['previous', 'first', 'mean']
        :param reference:
            *  *'previous'*: Aligns each frame (image) to its previous frame in the
               stack
            *  *'first:'* Aligns each frame (image) to the first frame in the stack -
               if n_frames is > 1, then each frame is  aligned to the
               mean of the first n_frames of the stack
            *  *'mean'*: Aligns each frame (image) to the average of all images in the
               stack

        :type n_frames: int, optional
        :param n_frames: If reference is 'first', then this parameter specifies the
                         number of frames from the beginning of the stack that should
                         be averaged to yield the reference image.

        :type axis: int, optional
        :param axis: The axis of the time dimension

        :type moving_average: int, optional
        :param moving_average:
            If moving_average is greater than 1, a moving average of the stack is first
            created (using a subset size of moving_average) before registration

        :type verbose: bool, optional
        :param verbose:
            Specifies whether a progressbar should be shown using tqdm.

        :type progress_callback: function, optional
        :param progress_callback:
            A function that is called after every iteration. This function should accept
            the keyword arguments current_iteration:int and end_iteration:int.
            It is always called from the thread that called this function,
            also when ``n_jobs`` > 1.

        :type suppress_axis_warning: bool, optional
        :param suppress_axis_warning:
            pystackreg will automatically try to determine the time axis and raise a
            warning when the detected time axis is not equal to the supplied axis.
            Set this option to True to suppress this warning.

        :type n_jobs: int, optional
        :param n_jobs:
            Number of threads used to register frames in parallel. 1 (default)
            registers sequentially in the calling thread; -1 uses all
            available CPUs. The result does not depend on ``n_jobs``.

        :rtype:  ndarray(img.shape[axis], 3, 3)
        :return: The transformation matrix for each image in the stack
        """
        stack, axis = self._check_stack(img, axis)
        n_workers = resolve_n_jobs(n_jobs)

        if reference not in self._valid_references:
            raise ValueError(f'Unknown reference "{reference}"')

        if self._transformation == self.BILINEAR and reference == "previous":
            raise ValueError(
                "Bilinear stack transformation not supported with reference "
                '== "previous" as a combination of bilinear transformations does not '
                "generally result in a bilinear transformation. Use another reference "
                "or manually register/transform images to their previous image."
            )

        n = stack.shape[axis]
        n_frames = _check_count(n_frames, "n_frames")
        # 0 is accepted and means "no moving average", like 1 (as in 0.2.x)
        moving_average = max(1, _check_count(moving_average, "moving_average", 0))
        if reference == "first" and n_frames > n:
            raise ValueError(
                f"n_frames ({n_frames}) exceeds the number of frames in the stack ({n})"
            )

        if not suppress_axis_warning:
            lowest_var_axis = self._detect_time_axis(stack)
            if axis != lowest_var_axis:
                warnings.warn(
                    f"Detected axis {lowest_var_axis} as the possible time axis for "
                    f"the stack due to its low variability, but axis {axis} was "
                    "supplied for registration. Are you sure you supplied the "
                    "correct axis?",
                    stacklevel=2,
                )

        frames = FrameSource(stack, axis, moving_average)

        # Frame 0 keeps the identity unless it is itself registered to a
        # reference image it is not identical to.
        first_index = 1
        ref: np.ndarray | None = None
        if reference == "first":
            ref = frames.mean_of_first(n_frames)
            if moving_average > 1:
                first_index = 0
        elif reference == "mean":
            ref = frames.mean_of_all()
            first_index = 0

        tmatdim = 4 if self._transformation == self.BILINEAR else 3
        tmats = np.repeat(
            np.identity(tmatdim).reshape((1, tmatdim, tmatdim)), n, axis=0
        ).astype(np.double)

        # Every task registers one frame to its reference and is independent
        # of all other tasks, also for reference="previous": there, the
        # pairwise transformations are chained after all tasks have finished.
        last: dict[str, Any] = {}

        def register_frame(i: int) -> None:
            reference_image = frames.frame(i - 1) if ref is None else ref
            m, refpts, movpts = self._register_pair(reference_image, frames.frame(i))
            tmats[i] = self._matrix_short_to_long(m)
            if i == n - 1:
                last["result"] = (m, refpts, movpts)

        indices = range(first_index, n)
        with Progress(len(indices), verbose, progress_callback, "Registering") as p:
            run_tasks(register_frame, indices, n_workers, p.step)

        if reference == "previous":
            for i in range(1, n):
                tmats[i] = np.matmul(tmats[i], tmats[i - 1])

        if "result" in last:
            self._m, self._refpts, self._movpts = last["result"]
            self._is_registered = True
        self._tmats = tmats

        return self._tmats

    def transform_stack(
        self,
        img: npt.ArrayLike,
        axis: int = 0,
        tmats: npt.ArrayLike | None = None,
        *,
        n_jobs: int | None = 1,
        out: np.ndarray | None = None,
        verbose: bool = False,
        progress_callback: Callable[..., Any] | None = None,
    ) -> np.ndarray:
        """
        Transform a stack after registration.

        :type img: array_like(Ni..., Nj..., Nk...)
        :param img: The 3D stack of images that should be aligned

        :type axis: int, optional
        :param axis: The axis of the time dimension

        :type tmats: ndarray(img.shape[axis], 3, 3), optional
        :param tmats: The transformation matrix for each image in the stack

        :type n_jobs: int, optional
        :param n_jobs:
            Number of threads used to transform frames in parallel. 1 (default)
            transforms sequentially in the calling thread; -1 uses all
            available CPUs. The result does not depend on ``n_jobs``.

        :type out: ndarray, optional
        :param out:
            Array of the same shape as ``img`` that receives the result. By
            default a new float64 array is allocated, which for integer input
            needs several times the memory of the stack. Pass an array of a
            smaller floating-point dtype, an integer dtype (values are rounded
            and clipped to the range of the dtype), a ``numpy.memmap``, or
            ``img`` itself to transform the stack in place.

        :type verbose: bool, optional
        :param verbose:
            Specifies whether a progressbar should be shown using tqdm.

        :type progress_callback: function, optional
        :param progress_callback:
            A function that is called after every transformed frame with the
            keyword arguments current_iteration:int and end_iteration:int.

        :rtype:  ndarray(Ni..., Nj..., Nk...)
        :return: The transformed stack (``out`` if it was given, otherwise a
                 new float64 array)
        """
        stack, axis = self._check_stack(img, axis)
        n_workers = resolve_n_jobs(n_jobs)

        if tmats is None:
            if self._tmats is None:
                raise RuntimeError(
                    "No transformation matrices given. Please register first or pass "
                    "transformation matrices explicitly."
                )
            matrices = self._tmats
        else:
            matrices = np.asarray(tmats)

        n = stack.shape[axis]
        tmatdim = 4 if self._transformation == self.BILINEAR else 3
        if matrices.ndim != 3 or matrices.shape[1:] != (tmatdim, tmatdim):
            raise ValueError(
                f"tmats must have shape (N, {tmatdim}, {tmatdim}), got {matrices.shape}"
            )
        if matrices.shape[0] != n:
            raise ValueError(
                "Number of saved transformation matrices does not match stack length"
            )

        if out is None:
            result = np.empty(stack.shape, dtype=np.float64)
        else:
            result = check_output_array(out, stack.shape)

        def transform_frame(i: int) -> None:
            frame = self.transform(simple_slice(stack, i, axis), matrices[i])
            store_frame(result, i, axis, frame)

        with Progress(n, verbose, progress_callback, "Transforming") as p:
            run_tasks(transform_frame, range(n), n_workers, p.step)

        return result

    def register_transform_stack(
        self,
        img: npt.ArrayLike,
        reference: str = "previous",
        n_frames: int = 1,
        axis: int = 0,
        moving_average: int = 1,
        verbose: bool = False,
        progress_callback: Callable[..., Any] | None = None,
        *,
        suppress_axis_warning: bool = False,
        n_jobs: int | None = 1,
        out: np.ndarray | None = None,
    ) -> np.ndarray:
        """
        Register and transform stack of images (movie).

        :type img: array_like(Ni..., Nj..., Nk...)
        :param img: The 3D stack of images that should be aligned

        :type reference: string, one of ['previous', 'first', 'mean']
        :param reference:
            *  *'previous'*: Aligns each frame (image) to its previous frame in the
               stack
            *  *'first:'* Aligns each frame (image) to the first frame in the stack -
               if n_frames is > 1, then each frame is aligned to the
               mean of the first n_frames of the stack
            *  *'mean'*: Aligns each frame (image) to the average of all images in the
               stack

        :type n_frames: int, optional
        :param n_frames: If reference is 'first', then this parameter specifies the
                         number of frames from the beginning of the stack that
                         should be averaged to yield the reference image.

        :type axis: int, optional
        :param axis: The axis of the time dimension

        :type moving_average: int, optional
        :param moving_average: If moving_average is greater than 1, a moving average of
                               the stack is first created (using a subset size of
                               moving_average) before registration

        :type verbose: bool, optional
        :param verbose:
            Specifies whether a progressbar should be shown using tqdm (one
            bar for the registration and one for the transformation).

        :type progress_callback: function, optional
        :param progress_callback:
            A function that is called after every iteration of the
            registration. This function should accept
            the keyword arguments current_iteration:int and end_iteration:int.

        :type suppress_axis_warning: bool, optional
        :param suppress_axis_warning:
            Set this option to True to suppress the warning that is raised
            when the supplied axis does not look like the time axis.

        :type n_jobs: int, optional
        :param n_jobs:
            Number of threads used for registration and transformation. 1
            (default) is sequential; -1 uses all available CPUs.

        :type out: ndarray, optional
        :param out:
            Array that receives the transformed stack, see
            :meth:`transform_stack`.

        :rtype:  ndarray(Ni..., Nj..., Nk...)
        :return: The transformed stack
        """
        stack = np.asarray(img)
        self.register_stack(
            stack,
            reference,
            n_frames,
            axis,
            moving_average,
            verbose,
            progress_callback,
            suppress_axis_warning,
            n_jobs=n_jobs,
        )

        return self.transform_stack(
            stack, axis, n_jobs=n_jobs, out=out, verbose=verbose
        )

    @staticmethod
    def _check_stack(img: npt.ArrayLike, axis: int) -> tuple[np.ndarray, int]:
        """Validate a stack argument and normalise ``axis`` to 0, 1 or 2."""
        stack = np.asarray(img)
        if stack.ndim != 3:
            raise ValueError(
                f"Stack must have three dimensions, got {stack.ndim} dimension(s)"
            )
        try:
            axis = operator.index(axis)
        except TypeError:
            raise TypeError(f"axis must be an integer, got {axis!r}") from None
        if not -3 <= axis < 3:
            raise ValueError(f"axis {axis} is out of bounds for a 3-D stack")
        axis %= 3
        if stack.shape[axis] == 0:
            raise ValueError(f"Stack has no frames along axis {axis}")
        return stack, axis


def _check_count(value: int, name: str, minimum: int = 1) -> int:
    """Validate an integer argument that must not be smaller than ``minimum``."""
    try:
        value = operator.index(value)
    except TypeError:
        raise TypeError(f"{name} must be an integer, got {value!r}") from None
    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}, got {value}")
    return value


def _bilinear_map(
    mat: np.ndarray, x: np.ndarray, y: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate a bilinear transformation (4x4 matrix acting on x, y, xy, 1)."""
    a, b = mat[0], mat[1]
    return (
        a[0] * x + a[1] * y + a[2] * x * y + a[3],
        b[0] * x + b[1] * y + b[2] * x * y + b[3],
    )


def _bilinear_solve(
    mat: np.ndarray, x: np.ndarray, y: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Invert a bilinear transformation for the given target coordinates.

    A bilinear map has no matrix inverse. The solution that is continuously
    connected to the affine case is found by Newton iteration started at the
    target coordinates.
    """
    a, b = mat[0], mat[1]
    u = np.array(x, dtype=np.float64, copy=True)
    v = np.array(y, dtype=np.float64, copy=True)
    for _ in range(100):
        fu = a[0] * u + a[1] * v + a[2] * u * v + a[3] - x
        fv = b[0] * u + b[1] * v + b[2] * u * v + b[3] - y
        j11, j12 = a[0] + a[2] * v, a[1] + a[2] * u
        j21, j22 = b[0] + b[2] * v, b[1] + b[2] * u
        det = j11 * j22 - j12 * j21
        with np.errstate(divide="ignore", invalid="ignore"):
            du = (j22 * fu - j12 * fv) / det
            dv = (j11 * fv - j21 * fu) / det
        u, v = u - du, v - dv
        step = np.max(np.hypot(du, dv)) if du.size else 0.0
        if not np.isfinite(step):
            break
        if step < 1e-10:
            return u, v
    raise RuntimeError(
        "Could not invert the bilinear transformation for the given points "
        "(the transformation is singular or too strongly distorted)"
    )
