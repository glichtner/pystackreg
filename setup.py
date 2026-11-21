# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2018-2026 Gregor Lichtner and PyStackReg contributors
"""Build script for the native ``pystackreg.turboreg`` extension.

All project metadata lives in ``pyproject.toml``; this file only describes the
C++ extension module, which setuptools cannot yet declare there in a stable way.
"""

import sys

import numpy
from setuptools import Extension, setup

if sys.platform == "win32":
    # MSVC: C++14 is the default (and minimum) language level; the default
    # floating-point model (/fp:precise) does not contract multiply-add.
    extra_compile_args = []
else:
    extra_compile_args = [
        "-std=c++11",
        # Forbid fusing a*b+c into a single FMA instruction. Contraction is on
        # by default on some platforms (e.g. arm64) and would make results
        # differ in the last bits between architectures.
        "-ffp-contract=off",
    ]

setup(
    ext_modules=[
        Extension(
            "pystackreg.turboreg",
            sources=[
                "src/pymain.cpp",
                "src/TurboReg.cpp",
                "src/TurboRegMask.cpp",
                "src/TurboRegImage.cpp",
                "src/TurboRegTransform.cpp",
                "src/TurboRegPointHandler.cpp",
            ],
            depends=[
                "inc/TurboReg.h",
                "inc/TurboRegImage.h",
                "inc/TurboRegMask.h",
                "inc/TurboRegPointHandler.h",
                "inc/TurboRegTransform.h",
                "inc/matrix.h",
            ],
            include_dirs=["inc", numpy.get_include()],
            extra_compile_args=extra_compile_args,
            language="c++",
        )
    ],
)
