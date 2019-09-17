#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "libaxb.h"
#include "libaxb/backend.h"
#include <assert.h>

#include <cuda_runtime.h>


static axbStatus_t cuda_malloc(void **ptr, size_t size_in_bytes, void *aux_data)
{
  (void)aux_data;
  cudaError_t err = cudaMalloc(ptr, size_in_bytes);
  //printf("Allocated device ptr %p\n", *ptr);
  assert(err == 0 && "cudaMalloc failed!");
  return err;
}

static axbStatus_t cuda_free(void *ptr_to_free, void *aux_data)
{
  (void)aux_data;
  cudaFree(ptr_to_free);
  return 0;
}


static axbStatus_t cuda_copyin(void *src, axbDataType_t src_type, void *dest, axbDataType_t dest_type, size_t n, void *aux_data)
{
  if (src_type != AXB_REAL_DOUBLE || dest_type != AXB_REAL_DOUBLE) return 17590; // not yet supported

  (void)aux_data;
  return cudaMemcpy(dest, src, sizeof(double) * n, cudaMemcpyHostToDevice);
}

static axbStatus_t cuda_copyout(void *src, axbDataType_t src_type, void *dest, axbDataType_t dest_type, size_t n, void *aux_data)
{
  if (src_type != AXB_REAL_DOUBLE || dest_type != AXB_REAL_DOUBLE) return 17590; // not yet supported

  (void)aux_data;
  return cudaMemcpy(dest, src, sizeof(double) * n, cudaMemcpyDeviceToHost);
}

extern "C" axbStatus_t axbMemBackendRegister_CUDA(struct axbHandle_s *handle)
{
  struct axbMemBackend_s *cuda_backend;
  axbMemBackendCreate(&cuda_backend);

  // populate host_backend:
  axbMemBackendSetName(cuda_backend, "CUDA");
  axbMemBackendSetMalloc(cuda_backend, cuda_malloc);
  axbMemBackendSetFree(cuda_backend, cuda_free);

  axbMemBackendSetCopyIn(cuda_backend, cuda_copyin);
  axbMemBackendSetCopyOut(cuda_backend, cuda_copyout);

  // push into enclosing context identified by handle:
  axbMemBackendRegister(handle, cuda_backend);
  return 0;
}

