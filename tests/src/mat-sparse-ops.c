
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libaxb.h"

#define TEST_EPSILON   1e-13
#include "common.h"



axbStatus_t test(axbMat_t A, axbMat_t B, axbVec_t x0, axbVec_t x1)
{
  // TODO: fill
  (void)A;
  (void)B;
  (void)x0;
  (void)x1;

  return 0;
}



axbStatus_t setup(axbHandle_t axb_handle, axbMemBackend_t mem, axbOpBackend_t ops)
{
  axbVec_t x0, x1;
  size_t n = 1297;

  // create vectors:
#define AXB_INIT_VEC(VECX) do { \
  AXB_ERR_CHECK(axbVecCreateBegin(axb_handle, &VECX)); \
  AXB_ERR_CHECK(axbVecSetSize(VECX, n)); \
  AXB_ERR_CHECK(axbVecSetMemBackend(VECX, mem)); \
  AXB_ERR_CHECK(axbVecSetOpBackend(VECX, ops)); \
  AXB_ERR_CHECK(axbVecCreateEnd(VECX)); \
} while (0)

  AXB_INIT_VEC(x0);
  AXB_INIT_VEC(x1);

  // create matrices:
  axbMat_t A;
  AXB_ERR_CHECK(axbMatCreateBegin(axb_handle, &A));
  AXB_ERR_CHECK(axbMatSetSizes(A, n, n));
  AXB_ERR_CHECK(axbMatSetMemBackend(A, mem));
  AXB_ERR_CHECK(axbMatSetOpBackend(A, ops));
  AXB_ERR_CHECK(axbMatSetStorageType(A, AXB_STORAGE_COMPRESSED_CSR));
  AXB_ERR_CHECK(axbMatCreateEnd(A));

  axbMat_t B;
  AXB_ERR_CHECK(axbMatCreateBegin(axb_handle, &B));
  AXB_ERR_CHECK(axbMatSetSizes(B, n, n));
  AXB_ERR_CHECK(axbMatSetMemBackend(B, mem));
  AXB_ERR_CHECK(axbMatSetOpBackend(B, ops));
  AXB_ERR_CHECK(axbMatSetStorageType(B, AXB_STORAGE_COMPRESSED_CSR));
  AXB_ERR_CHECK(axbMatCreateEnd(B));

  // Run test
  test(A, B, x0, x1);

  AXB_ERR_CHECK(axbMatDestroy(A));
  AXB_ERR_CHECK(axbMatDestroy(B));

  AXB_ERR_CHECK(axbVecDestroy(x0));
  AXB_ERR_CHECK(axbVecDestroy(x1));

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
