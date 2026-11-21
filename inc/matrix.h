/*
 * This file is part of PyStackReg.
 *
 * Copyright (C) 2018-2026 Gregor Lichtner and PyStackReg contributors
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version. It is distributed WITHOUT ANY
 * WARRANTY; see the GNU General Public License for more details
 * (file LICENSE or <https://www.gnu.org/licenses/>).
 */

/* Minimal row-major matrix used in place of the Java double[][] arrays. */

#include <vector>

#ifndef MATRIX_H_
#define MATRIX_H_

template <typename T> class matrix {
  std::vector<T> m;

public:
  matrix() {}

  matrix(unsigned int rows, unsigned int cols) { resize(rows, cols); }

  void resize(unsigned int rows, unsigned int cols) {
    this->rows = rows;
    this->cols = cols;
    m.resize(rows * cols);
  }

  T &operator()(unsigned int const i, unsigned int const j) {
    return m[i * cols + j];
  }

  T *getRowPtr(unsigned int const i) { return &m[0] + i * cols; }

  T *begin() { return &m[0]; }

  T *end() { return &m[0] + rows * cols; }

  T *getPtr(unsigned int const i, unsigned int const j) {
    return &m[0] + i * cols + j;
  }

  unsigned int nrows() { return rows; }

  unsigned int ncols() { return cols; }

private:
  unsigned int cols;
  unsigned int rows;
};
#endif
