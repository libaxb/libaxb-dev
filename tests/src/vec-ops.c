
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libaxb.h"

#define TEST_EPSILON   1e-13
#include "common.h"



axbStatus_t test(axbVec_t x0, axbVec_t x1, axbVec_t x2, axbVec_t x3, axbScalar_t alpha, axbScalar_t beta, axbScalar_t gamma,
                 double *y0, double *y1, double *y2, double *y3, double a, double b, double c)
{
  size_t n;
  AXB_ERR_CHECK(axbVecGetSize(x1, &n));
  double *temp = malloc(sizeof(double) * n);
  double t = 0;
  double s = 0;



  //
  // in-place operations
  //

  // AXB_OP_VEC_SET
  AXB_ERR_CHECK(axbVecSet(x0, alpha));
  for (size_t i=0; i<n; ++i) y0[i] = a;
  check_equal(x0, y0, temp, "axbVecSet");

  // AXB_OP_VEC_SQRTABS
  AXB_ERR_CHECK(axbVecSqrtAbs(x0));
  for (size_t i=0; i<n; ++i) y0[i] = sqrt(fabs(y0[i]));
  check_equal(x0, y0, temp, "axbVecSqrtAbs");

  // AXB_OP_VEC_ZERO
  AXB_ERR_CHECK(axbVecZero(x0));
  for (size_t i=0; i<n; ++i) y0[i] = 0;
  check_equal(x0, y0, temp, "axbVecZero");

  // AXB_OP_VEC_SCALE
  AXB_ERR_CHECK(axbVecScale(x1, alpha));
  for (size_t i=0; i<n; ++i) y1[i] *= a;
  check_equal(x1, y1, temp, "axbVecScale");



  //
  // reduction operations
  //

  //AXB_OP_VEC_SUM
  AXB_ERR_CHECK(axbVecSum(x1, alpha));
  t=0;
  for (size_t i=0; i<n; ++i) t += y1[i];
  check_equal_scalar(alpha, t, "axbVecSum");


  //AXB_OP_VEC_DOT
  AXB_ERR_CHECK(axbVecDot(x1, x2, alpha));
  t=0;
  for (size_t i=0; i<n; ++i) t += y1[i] * y2[i];
  check_equal_scalar(alpha, t, "axbVecDot");

  //AXB_OP_VEC_TDOT
  AXB_ERR_CHECK(axbVecTDot(x1, x2, alpha));
  t=0;
  for (size_t i=0; i<n; ++i) t += y1[i] * y2[i];
  check_equal_scalar(alpha, t, "axbVecTDot");

  //AXB_OP_VEC_MDOT
  // skipping: tested separately

  //AXB_OP_VEC_NORM1
  AXB_ERR_CHECK(axbVecNorm1(x1, alpha));
  t=0;
  for (size_t i=0; i<n; ++i) t += fabs(y1[i]);
  check_equal_scalar(alpha, t, "axbVecNorm1");

  //AXB_OP_VEC_NORM2
  AXB_ERR_CHECK(axbVecNorm2(x1, alpha));
  t=0;
  for (size_t i=0; i<n; ++i) t += y1[i] * y1[i];
  check_equal_scalar(alpha, sqrt(t), "axbVecNorm2");

  //AXB_OP_VEC_NORMINF
  AXB_ERR_CHECK(axbVecNormInf(x1, alpha));
  t=0;
  for (size_t i=0; i<n; ++i) t = fabs(y1[i]) > t ? fabs(y1[i]) : t;
  check_equal_scalar(alpha, t, "axbVecNormInf");

  //AXB_OP_VEC_DOTNORM2
  AXB_ERR_CHECK(axbVecDotNorm2(x1, x2, alpha, beta));
  t=0;
  s=0;
  for (size_t i=0; i<n; ++i) { s += y1[i] * y2[i]; t += y2[i] * y2[i]; }
  check_equal_scalar(alpha, s, "axbVecDotNorm2, alpha");
  check_equal_scalar(beta, sqrt(t), "axbVecDotNorm2, beta");

  //AXB_OP_VEC_MAX
  size_t idx0 = 0;
  AXB_ERR_CHECK(axbVecMax(x1, &idx0, alpha));
  t=y1[0];
  size_t idx1 = 0;
  for (size_t i=1; i<n; ++i) {
    if (y1[i] > t) { idx1 = i; t = y1[i]; }
  }
  AXB_ERR_CHECK(axbScalarGetValue(alpha, &s, AXB_REAL_DOUBLE));
  if (idx0 != idx1) { fprintf(stderr, "ERROR: index mismatch for axbVecMax: %ld (%g) vs. %ld (%g)\n", idx0, s, idx1, t); exit(EXIT_FAILURE); }
  check_equal_scalar(alpha, t, "axbVecMax");

  //AXB_OP_VEC_MIN
  idx0 = 0;
  AXB_ERR_CHECK(axbVecMin(x1, &idx0, alpha));
  t=y1[0];
  idx1 = 0;
  for (size_t i=1; i<n; ++i) {
    if (y1[i] < t) { idx1 = i; t = y1[i]; }
  }
  AXB_ERR_CHECK(axbScalarGetValue(alpha, &s, AXB_REAL_DOUBLE));
  if (idx0 != idx1) { fprintf(stderr, "ERROR: index mismatch for axbVecMin: %ld (%g) vs. %ld (%g)\n", idx0, s, idx1, t); exit(EXIT_FAILURE); }
  check_equal_scalar(alpha, t, "axbVecMin");



  //
  // vector-vector operations
  //

  //AXB_OP_VEC_COPY
  AXB_ERR_CHECK(axbVecCopy(x1, x0));
  for (size_t i=0; i<n; ++i) y0[i] = y1[i];
  check_equal(x0, y0, temp, "axbVecCopy, x0");
  check_equal(x1, y1, temp, "axbVecCopy, x1");

  //AXB_OP_VEC_SWAP
  AXB_ERR_CHECK(axbVecSwap(x0, x2));
  for (size_t i=0; i<n; ++i) { t = y0[i]; y0[i] = y2[i]; y2[i] = t; }
  check_equal(x0, y0, temp, "axbVecSwap, x0");
  check_equal(x2, y2, temp, "axbVecSwap, x2");

  AXB_ERR_CHECK(axbScalarGetValue(alpha, &a, AXB_REAL_DOUBLE));  // reset
  AXB_ERR_CHECK(axbScalarGetValue(beta,  &b, AXB_REAL_DOUBLE));  // reset
  AXB_ERR_CHECK(axbScalarGetValue(gamma, &c, AXB_REAL_DOUBLE));  // reset

  //AXB_OP_VEC_AXPY
  AXB_ERR_CHECK(axbVecAXPY(x0, alpha, x1));
  for (size_t i=0; i<n; ++i) y0[i] += a * y1[i];
  check_equal(x0, y0, temp, "axbVecAXPY");

  //AXB_OP_VEC_AYPX
  AXB_ERR_CHECK(axbVecAYPX(x0, alpha, x1));
  for (size_t i=0; i<n; ++i) y0[i] = a * y0[i] + y1[i];
  check_equal(x0, y0, temp, "axbVecAYPX");

  //AXB_OP_VEC_AXPBYPCZ
  AXB_ERR_CHECK(axbVecAXPBYPCZ(x2, alpha, beta, gamma, x0, x1));
  for (size_t i=0; i<n; ++i) y2[i] = a * y0[i] + b * y1[i] + c * y2[i];
  check_equal(x2, y2, temp, "axbVecAXPBYPCZ");

  //AXB_OP_VEC_WAXPY
  AXB_ERR_CHECK(axbVecWAXPY(x2, alpha, x0, x1));
  for (size_t i=0; i<n; ++i) y2[i] = a * y0[i] + y1[i];
  check_equal(x2, y2, temp, "axbVecWAXPY");

  //AXB_OP_VEC_MAXPY
  axbVec_t vecs[3] = {x1, x2, x3};
  axbScalar_t scalars[3] = {alpha, beta, gamma};
  AXB_ERR_CHECK(axbVecMAXPY(x0, 3, scalars, vecs));
  for (size_t i=0; i<n; ++i) y0[i] = a * y1[i] + b * y2[i] + c * y3[i];
  check_equal(x0, y0, temp, "axbVecMAXPY");

  //AXB_OP_VEC_POINTWISEMULT
  AXB_ERR_CHECK(axbVecPointwiseMult(x0, x1, x2));
  for (size_t i=0; i<n; ++i) y0[i] = y1[i] * y2[i];
  check_equal(x0, y0, temp, "axbVecPointwiseMult");

  //AXB_OP_VEC_POINTWISEDIV
  AXB_ERR_CHECK(axbVecPointwiseDivide(x0, x1, x2));
  for (size_t i=0; i<n; ++i) y0[i] = y1[i] / y2[i];
  check_equal(x0, y0, temp, "axbVecPointwiseDivide");

  free(temp);
  return 0;
}


axbStatus_t setup(axbHandle_t axb_handle, axbMemBackend_t mem, axbOpBackend_t ops)
{
  axbVec_t x0, x1, x2, x3;
  double *y0, *y1, *y2, *y3;
  size_t n = 1297;

  // create vectors:
#define AXB_INIT_VECS(VECX, VECY) do { \
  AXB_ERR_CHECK(axbVecCreateBegin(axb_handle, &VECX)); \
  AXB_ERR_CHECK(axbVecSetSize(VECX, n)); \
  AXB_ERR_CHECK(axbVecSetMemBackend(VECX, mem)); \
  AXB_ERR_CHECK(axbVecSetOpBackend(VECX, ops)); \
  AXB_ERR_CHECK(axbVecCreateEnd(VECX)); \
  VECY = malloc(sizeof(double) * n); \
  initVecRandom(VECX, VECY); } while (0)

  AXB_INIT_VECS(x0, y0);
  AXB_INIT_VECS(x1, y1);
  AXB_INIT_VECS(x2, y2);
  AXB_INIT_VECS(x3, y3);

  // create scalars:
  axbScalar_t alpha;
  double a = 3.1415;
  AXB_ERR_CHECK(axbScalarCreate(axb_handle, &alpha, &a, AXB_REAL_DOUBLE, mem));

  axbScalar_t beta;
  double b = 2.717;
  AXB_ERR_CHECK(axbScalarCreate(axb_handle, &beta, &b, AXB_REAL_DOUBLE, mem));

  axbScalar_t gamma;
  double c = 1.4142;
  AXB_ERR_CHECK(axbScalarCreate(axb_handle, &gamma, &c, AXB_REAL_DOUBLE, mem));

  // Run test
  test(x0, x1, x2, x3, alpha, beta, gamma,
       y0, y1, y2, y3, a,     b,    c);

  AXB_ERR_CHECK(axbVecDestroy(x0)); free(y0);
  AXB_ERR_CHECK(axbVecDestroy(x1)); free(y1);
  AXB_ERR_CHECK(axbVecDestroy(x2)); free(y2);
  AXB_ERR_CHECK(axbVecDestroy(x3)); free(y3);
  AXB_ERR_CHECK(axbScalarDestroy(alpha));
  AXB_ERR_CHECK(axbScalarDestroy(beta));
  AXB_ERR_CHECK(axbScalarDestroy(gamma));

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
