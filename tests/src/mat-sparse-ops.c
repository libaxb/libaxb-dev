
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libaxb.h"

#define TEST_EPSILON   1e-13
#include "common.h"



axbStatus_t test(axbMat_t A, axbMat_t B, axbVec_t x0, axbVec_t x1)
{
  size_t rows, cols, nonzeros;
  AXB_ERR_CHECK(axbMatGetSizes(A, &rows, &cols));
  AXB_ERR_CHECK(axbMatGetNonzerosSize(A, &nonzeros));
  int *A_row_markers = malloc(sizeof(int) * (rows+1));
  int *A_col_indices = malloc(sizeof(int) * nonzeros);
  double *A_values   = malloc(sizeof(double) * nonzeros);

  double *temp = malloc(sizeof(double) * (rows+1));
  double *x_host = malloc(sizeof(double) * (rows+1));
  double *y_host = malloc(sizeof(double) * (rows+1));

  AXB_ERR_CHECK(axbMatGetValuesCSR(A, A_row_markers, AXB_INT_32, A_col_indices, AXB_INT_32, A_values, AXB_REAL_DOUBLE));
  AXB_ERR_CHECK(axbVecGetValues(x0, x_host, AXB_REAL_DOUBLE));

  // AXB_OP_MAT_VEC
  AXB_ERR_CHECK(axbMatVec(A, x0, x1));
  for (size_t i=0; i<rows; ++i) {  // reference calculation
    double val = 0;
    for (int j=A_row_markers[i]; j<A_row_markers[i+1]; ++j) {
      val += A_values[j] * x_host[A_col_indices[j]];
    }
    y_host[i] = val;
  }
  check_equal(x1, y_host, temp, "axbMatVec");

  // AXB_OP_MAT_TVEC
  AXB_ERR_CHECK(axbMatTVec(A, x0, x1));
  for (size_t i=0; i<cols; ++i) y_host[i] = 0;
  for (size_t i=0; i<rows; ++i) {  // reference calculation
    for (int j=A_row_markers[i]; j<A_row_markers[i+1]; ++j) {
      y_host[A_col_indices[j]] += A_values[j] * x_host[i];
    }
  }
  check_equal(x1, y_host, temp, "axbMatTVec");

  // Skipping matrix-matrix multiplication test

  // AXB_OP_MAT_TVEC
  AXB_ERR_CHECK(axbMatTVec(A, x0, x1));
  for (size_t i=0; i<cols; ++i) y_host[i] = 0;
  for (size_t i=0; i<rows; ++i) {  // reference calculation
    for (int j=A_row_markers[i]; j<A_row_markers[i+1]; ++j) {
      y_host[A_col_indices[j]] += A_values[j] * x_host[i];
    }
  }
  check_equal(x1, y_host, temp, "axbMatTVec");


  // AXB_OP_MAT_TRANS
  axbMat_t C, D;
  AXB_ERR_CHECK(axbMatTrans(A, &C));
  AXB_ERR_CHECK(axbMatVec(C, x0, x1));
  check_equal(x1, y_host, temp, "axbMatTrans (via axbMatVec)");
  AXB_ERR_CHECK(axbMatTrans(C, &D));
  check_equal(x1, y_host, temp, "axbMatTrans (double-transposition)");


  //
  // clean up:
  //
  AXB_ERR_CHECK(axbMatDestroy(D));
  AXB_ERR_CHECK(axbMatDestroy(C));

  free(A_row_markers);
  free(A_col_indices);
  free(A_values);

  free(temp);
  free(x_host);
  free(y_host);

  return 0;
}



axbStatus_t setup(axbHandle_t axb_handle, axbMemBackend_t mem, axbOpBackend_t ops)
{
  axbVec_t x0, x1;
  size_t n = 154534;
  double *temp = malloc(sizeof(double) * n);

  // create vectors:
#define AXB_INIT_VEC(VECX) do { \
  AXB_ERR_CHECK(axbVecCreateBegin(axb_handle, &VECX)); \
  AXB_ERR_CHECK(axbVecSetSize(VECX, n)); \
  AXB_ERR_CHECK(axbVecSetMemBackend(VECX, mem)); \
  AXB_ERR_CHECK(axbVecSetOpBackend(VECX, ops)); \
  AXB_ERR_CHECK(axbVecCreateEnd(VECX)); \
  initVecRandom(VECX, temp); \
} while (0)

  AXB_INIT_VEC(x0);
  AXB_INIT_VEC(x1);

  // create matrices:
  axbMat_t A;
  AXB_ERR_CHECK(axbMatCreateBegin(axb_handle, &A));
  AXB_ERR_CHECK(axbMatSetSizes(A, n, n));
  AXB_ERR_CHECK(axbMatSetMemBackend(A, mem));
  AXB_ERR_CHECK(axbMatSetOpBackend(A, ops));
  AXB_ERR_CHECK(axbMatSetStorageType(A, AXB_STORAGE_CSR));
  AXB_ERR_CHECK(axbMatCreateEnd(A));
  AXB_ERR_CHECK(initMatSparseRandom(A, 5));

  axbMat_t B;
  AXB_ERR_CHECK(axbMatCreateBegin(axb_handle, &B));
  AXB_ERR_CHECK(axbMatSetSizes(B, n, n));
  AXB_ERR_CHECK(axbMatSetMemBackend(B, mem));
  AXB_ERR_CHECK(axbMatSetOpBackend(B, ops));
  AXB_ERR_CHECK(axbMatSetStorageType(B, AXB_STORAGE_CSR));
  AXB_ERR_CHECK(axbMatCreateEnd(B));
  AXB_ERR_CHECK(initMatSparseRandom(B, 5));

  // Run test
  test(A, B, x0, x1);

  AXB_ERR_CHECK(axbMatDestroy(A));
  AXB_ERR_CHECK(axbMatDestroy(B));

  AXB_ERR_CHECK(axbVecDestroy(x0));
  AXB_ERR_CHECK(axbVecDestroy(x1));

  free(temp);

  return 0;
}



int main(int argc, char **argv)
{
  axbHandle_t axbHandle;
  const char *mem_backend_name = argv[1];
  const char *op_backend_name = argv[2];

  if (argc < 2) {
    fprintf(stderr, "Missing argument: MemoryBackend\n");
    return EXIT_FAILURE;
  }
  else if (argc == 2) { // mem backend equals op backend
    op_backend_name = argv[1];
  }

  printf("Initializing libaxb...\n");
  AXB_ERR_CHECK(axbInit(&axbHandle));

  axbMemBackend_t mem_backend;
  AXB_ERR_CHECK(axbMemBackendGetByName(axbHandle, &mem_backend, mem_backend_name));
  if (!mem_backend) {
    fprintf(stderr, "ERROR: Cannot find memory backend %s\n", mem_backend_name);
    return EXIT_FAILURE;
  }
  printf("Using memory backend %s\n", mem_backend_name);

  axbOpBackend_t op_backend;
  AXB_ERR_CHECK(axbOpBackendGetByName(axbHandle, &op_backend, op_backend_name));
  if (!op_backend) {
    fprintf(stderr, "ERROR: Cannot find memory backend %s\n", mem_backend_name);
    return EXIT_FAILURE;
  }
  printf("Using op backend %s\n", op_backend_name);

  setup(axbHandle, mem_backend, op_backend);


  printf("Shutting down libaxb...\n");
  AXB_ERR_CHECK(axbFinalize(axbHandle));

  printf("Completed successfully\n");

  return EXIT_SUCCESS;
}
