# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2018-2026 Gregor Lichtner and PyStackReg contributors
"""pystackreg: image registration with a Python/C++ port of TurboReg/StackReg.

The registration algorithm is the one of the ImageJ plugin TurboReg by
Philippe Thévenaz (EPFL), described in

    P. Thévenaz, U.E. Ruttimann, M. Unser, "A Pyramid Approach to Subpixel
    Registration Based on Intensity", IEEE Transactions on Image Processing,
    vol. 7, no. 1, pp. 27-41, January 1998.
"""

from pystackreg.pystackreg import StackReg
from pystackreg.version import __version__

__all__ = ["StackReg", "__version__"]
