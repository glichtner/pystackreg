/*
 * This file is part of PyStackReg.
 *
 * It is a C++ translation of parts of the TurboReg ImageJ plugin by
 * Philippe Thevenaz (Java source file TurboReg_.java, version "July 7,
 * 2011"), with modifications. See NOTICE.
 *
 * Copyright (C) 2011 Philippe Thevenaz
 *     Original implementation (TurboReg).
 * Copyright (C) 2018-2026 Gregor Lichtner and PyStackReg contributors
 *     C++ translation and subsequent modifications.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/*====================================================================
| Philippe Thevenaz
| EPFL/STI/IMT/LIB/BM.4.137
| Station 17
| CH-1015 Lausanne VD
| Switzerland
|
| phone (CET): +41(21)693.51.61
| fax: +41(21)693.37.01
| RFC-822: philippe.thevenaz@epfl.ch
| X-400: /C=ch/A=400net/P=switch/O=epfl/S=thevenaz/G=philippe/
| URL: http://bigwww.epfl.ch/
\===================================================================*/

/*====================================================================
| This work is based on the following paper:
|
| P. Thevenaz, U.E. Ruttimann, M. Unser
| A Pyramid Approach to Subpixel Registration Based on Intensity
| IEEE Transactions on Image Processing
| vol. 7, no. 1, pp. 27-41, January 1998.
|
| This paper is available on-line at
| http://bigwww.epfl.ch/publications/thevenaz9801.html
|
| Other relevant on-line publications are available at
| http://bigwww.epfl.ch/publications/
\===================================================================*/

/*====================================================================
| Historical notice. The translated source file (TurboReg, version
| July 7, 2011) carried the notice quoted below. It is kept for the
| record and does NOT state the terms under which this file is
| distributed: its author has since released TurboReg under the GNU
| General Public License, version 3 or (at your option) any later
| version (https://bigwww.epfl.ch/thevenaz/turboreg/, "Conditions of
| Use"), and current versions of TurboReg_.java no longer carry it.
|
|   Additional help available at http://bigwww.epfl.ch/thevenaz/turboreg/
|
|   You'll be free to use this software for research purposes, but you
|   should not redistribute it without our consent. In addition, we expect
|   you to include a citation or acknowledgment whenever you present or
|   publish results that are based on it.
\===================================================================*/

#include <cmath>

#include "TurboReg.h"
#include "TurboRegImage.h"
#include "matrix.h"

#ifndef TURBOREGPOINTHANDLER_H_
#define TURBOREGPOINTHANDLER_H_

class TurboRegPointHandler { /* class TurboRegPointHandler */
public:
  TurboRegPointHandler(TurboRegImage &img, Transformation transformation);

  TurboRegPointHandler(matrix<double> &precisionPoint);

  matrix<double> &getPoints();

  void setPoints(matrix<double> &precisionPoint);

  void setPointsByTransformation(int width, int height,
                                 Transformation transformation);

private:
  /*********************************************************************
   The golden ratio mathematical constant determines where to put the
  initial landmarks.
  ********************************************************************/
  static double GOLDEN_RATIO;

  matrix<double> precisionPoint; // = new double[NUM_POINTS * 2];
};

#endif
