#ifndef AXB_TESTS_COMMON_H
#define AXB_TESTS_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include "libaxb.h"

#ifndef TEST_EPSILON
#define TEST_EPSILON   1e-13
#endif

#define AXB_ERR_CHECK( arg )   do { int err_check_12388123 = arg; if (err_check_12388123) { fprintf(stderr, "Failure in line %d in file %s\n", __LINE__, __FILE__); return err_check_12388123; }  } while (0)

double get_rel_diff(double a, double b) {
  double fmax = fabs(a) > fabs(b) ? fabs(a) : fabs(b);

  double rel_diff = (a - b); // also deals with all-zero result correctly
  if (fabs(rel_diff) > 0) rel_diff /= fmax;
  return rel_diff;
}


axbStatus_t check_equal_scalar(axbScalar_t alpha, double a, const char *test_desc)
{
  double t;
  AXB_ERR_CHECK(axbScalarGetValue(alpha, &t, AXB_REAL_DOUBLE));

  double rel_diff = get_rel_diff(t, a);
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
    double rel_diff = get_rel_diff(y[i], temp[i]);
    if (rel_diff > TEST_EPSILON) {
      fprintf(stderr, "ERROR: Test %s failed with relative error %g at index %ld: reference %g vs. libaxb %g\n Aborting...\n", test_desc, rel_diff, i, y[i], temp[i]);
      exit(EXIT_FAILURE);
    }
  }
  printf("PASSED: Test %s\n", test_desc);
  return 0;
}

axbStatus_t check_equal_csr(axbMat_t A, axbMat_t B, const char *test_desc)
{
  //
  // Check 1: Dimensions
  //
  size_t A_rows, A_cols, A_nnz;
  AXB_ERR_CHECK(axbMatGetSizes(A, &A_rows, &A_cols));
  AXB_ERR_CHECK(axbMatGetNonzerosSize(A, &A_nnz));

  size_t B_rows, B_cols, B_nnz;
  AXB_ERR_CHECK(axbMatGetSizes(B, &B_rows, &B_cols));
  AXB_ERR_CHECK(axbMatGetNonzerosSize(B, &B_nnz));

  //
  // Check 2: Values
  //
  int *A_row_markers = malloc(sizeof(int)*(A_rows+1));
  int *A_col_indices = malloc(sizeof(int)*A_nnz);
  double *A_values   = malloc(sizeof(double)*A_nnz);
  AXB_ERR_CHECK(axbMatGetValuesCSR(A, A_row_markers, AXB_INT_32, A_col_indices, AXB_INT_32, A_values, AXB_REAL_DOUBLE));


  int *B_row_markers = malloc(sizeof(int)*(B_rows+1));
  int *B_col_indices = malloc(sizeof(int)*B_nnz);
  double *B_values   = malloc(sizeof(double)*B_nnz);
  AXB_ERR_CHECK(axbMatGetValuesCSR(B, B_row_markers, AXB_INT_32, B_col_indices, AXB_INT_32, B_values, AXB_REAL_DOUBLE));


  // test row array:
  for (size_t i=0; i<A_rows; ++i) {
    double rel_diff = get_rel_diff(A_row_markers[i], B_row_markers[i]);

    if (rel_diff > TEST_EPSILON) {
      fprintf(stderr, "ERROR: Test %s (row array) failed with relative error %g: reference %d vs. libaxb %d\n Aborting...\n", test_desc, rel_diff, A_row_markers[i], B_row_markers[i]);
      exit(EXIT_FAILURE);
    }
  }

  // test col indices and values arrays:
  for (size_t i=0; i<A_nnz; ++i) {
    double rel_diff = get_rel_diff(A_col_indices[i], B_col_indices[i]);

    if (rel_diff > TEST_EPSILON) {
      fprintf(stderr, "ERROR: Test %s (col array) failed with relative error %g: reference %d vs. libaxb %d\n Aborting...\n", test_desc, rel_diff, A_col_indices[i], B_col_indices[i]);
      exit(EXIT_FAILURE);
    }

    rel_diff = get_rel_diff(A_values[i], B_values[i]);

    if (rel_diff > TEST_EPSILON) {
      fprintf(stderr, "ERROR: Test %s (col array) failed with relative error %g: reference %g vs. libaxb %g\n Aborting...\n", test_desc, rel_diff, A_values[i], B_values[i]);
      exit(EXIT_FAILURE);
    }

  }


  printf("PASSED: Test %s\n", test_desc);

  free(A_row_markers);
  free(A_col_indices);
  free(A_values);
  free(B_row_markers);
  free(B_col_indices);
  free(B_values);
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

  printf("initMatSparseRandom(): Initializing sparse matrix with %ld rows and %ld nonzeros\n", rows, nonzeros);

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

  free(host_rows);
  free(host_cols);
  free(host_vals);
  return 0;
}


axbStatus_t printMat(axbMat_t A)
{
  const char *name;
  axbMatGetName(A, &name);
  printf("Matrix %s:\n", name);

  size_t A_rows, A_cols, A_nnz;
  AXB_ERR_CHECK(axbMatGetSizes(A, &A_rows, &A_cols));
  AXB_ERR_CHECK(axbMatGetNonzerosSize(A, &A_nnz));

  int *A_row_markers = malloc(sizeof(int)*(A_rows+1));
  int *A_col_indices = malloc(sizeof(int)*A_nnz);
  double *A_values   = malloc(sizeof(double)*A_nnz);
  AXB_ERR_CHECK(axbMatGetValuesCSR(A, A_row_markers, AXB_INT_32, A_col_indices, AXB_INT_32, A_values, AXB_REAL_DOUBLE));

  for (size_t i=0; i<A_rows; ++i) {
    printf("Row %ld: ", i);
    for (int j=A_row_markers[i]; j<A_row_markers[i+1]; ++j) {
      printf("(%d,%g) ", A_col_indices[j], A_values[j]);
    }
    printf("\n");
  }
  printf("-----\n");

  return 0;
}

#endif
