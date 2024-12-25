#include <iostream>

#include "Kokkos_Core.hpp"

#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_Driver.hpp"

int main(int argc, char *argv[]) {

  using real = float;

  int n = 64;
  int nnz = 3 * n - 2;

  real h = real(1) / n;

  Kokkos::View<real*> values("values", nnz);
  Kokkos::View<int*> rowPtr("row pointer", n + 1), colIdx("column indices", nnz);

  rowPtr(0) = 0;
  rowPtr(1) = 2;
  colIdx(0) = 0; colIdx(1) = 1;
  values(0) = real(2.0 / h); values(1) = -real(1.0 / h);
  for (int i = 1; i < n; ++i) {
    colIdx(rowPtr(i)) = i - 1;
    colIdx(rowPtr(i) + 1) = i;
    values(rowPtr(i)) = -real(1.0 / h);
    values(rowPtr(i) + 1) = real(2.0 / h);
    if (i + 1 < n) {
    rowPtr(i + 1) = rowPtr(i) + 3;
    colIdx(rowPtr(i) + 2) = i + 1;
    values(rowPtr(i) + 2) = -real(1.0 / h);
    }
    else {
      rowPtr(i + 1) = rowPtr(i) + 2;
    }
  }

  return 0;
}

