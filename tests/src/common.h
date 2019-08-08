#ifndef AXB_TESTS_COMMON_H
#define AXB_TESTS_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include "libaxb.h"

#ifndef TEST_EPSILON
#define TEST_EPSILON   1e-13
#endif

#define AXB_ERR_CHECK( arg )   do { int err_check_12388123 = arg; if (err_check_12388123) { fprintf(stderr, "Failure in line %d in file %s\n", __LINE__, __FILE__); return err_check_12388123; }  } while (0)

axbStatus_t check_equal_scalar(axbScalar_t alpha, double a, const char *test_desc)
{
  double t;
  AXB_ERR_CHECK(axbScalarGetValue(alpha, &t, AXB_REAL_DOUBLE));

  double fmax = fabs(t);
  if (fmax < fabs(a)) fmax = fabs(a);
  double rel_diff = (t - a); // also deals with all-zero result correctly
  if (fabs(rel_diff) > 0) rel_diff /= fmax;

  if (rel_diff > TEST_EPSILON) {
    fprintf(stderr, "ERROR: Test %s failed with relative error %g: reference %g vs. libaxb %g\n Aborting...\n", test_desc, rel_diff, a, t);
    exit(EXIT_FAILURE);
  }

  printf("PASSED: Test %s\n", test_desc);
  return 0;
}

axbStatus_t check_equal(axbVec_t x, double *y, double *temp, const char *test_desc)
{
  size_t n;
  AXB_ERR_CHECK(axbVecGetSize(x, &n));

  AXB_ERR_CHECK(axbVecGetValues(x, temp, AXB_REAL_DOUBLE));
  for (size_t i=0; i<n; ++i) {
    double fmax = fabs(temp[i]);
    if (fmax < fabs(y[i])) fmax = fabs(y[i]);
    double rel_diff = (temp[i] - y[i]); // also deals with all-zero result correctly
    if (fabs(rel_diff) > 0) rel_diff /= fmax;

    if (rel_diff > TEST_EPSILON) {
      fprintf(stderr, "ERROR: Test %s failed with relative error %g: reference %g vs. libaxb %g\n Aborting...\n", test_desc, rel_diff, y[i], temp[i]);
      exit(EXIT_FAILURE);
    }
  }
  printf("PASSED: Test %s\n", test_desc);
  return 0;
}


axbStatus_t initVecRandom(axbVec_t x, double *y)
{
  size_t n;
  AXB_ERR_CHECK(axbVecGetSize(x, &n));

  // populate host vector:
  for (size_t i=0; i<n; ++i) {
    double r = rand();
    r /= RAND_MAX;
    y[i] = -2.0 + r;  // stay away from zero to avoid singularity issues
  }

  // copy over to axbVec:
  AXB_ERR_CHECK(axbVecSetValues(x, y, AXB_REAL_DOUBLE));
  return 0;
}


axbStatus_t initMatSparseRandom(axbMat_t A, size_t nonzeros_per_row)
{
  size_t rows, cols;
  AXB_ERR_CHECK(axbMatGetSizes(A, &rows, &cols));
  size_t nonzeros = rows * nonzeros_per_row;

  int    *host_rows = malloc(sizeof(int)   *(rows+1));
  int    *host_cols = malloc(sizeof(int)   *(nonzeros));
  double *host_vals = malloc(sizeof(double)*(nonzeros));

  // populate
  host_rows[0] = 0;
  int nnz_index = 0;
  size_t col_inc = rows / nonzeros_per_row / 2;
  for (size_t i=0; i<rows; ++i) {
    for (size_t j=0; j<nonzeros_per_row; ++j) {
      host_cols[nnz_index] = i/2 + col_inc * j + rand() % col_inc;
      host_vals[nnz_index] = 1.0 + (rand() % 1023 - 500.0) / 300.0;
      ++nnz_index;
    }
    host_rows[i+1] = nnz_index;
  }

  // copy over to matrix:
  AXB_ERR_CHECK(axbMatSetValuesCSR(A, host_rows, AXB_INT_32, host_cols, AXB_INT_32, host_vals, AXB_REAL_DOUBLE, nonzeros));
  return 0;
}

#endif
