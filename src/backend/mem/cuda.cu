#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "libaxb.h"
#include "libaxb/backend.h"

#include <cuda_runtime.h>


static void * cuda_malloc(size_t size_in_bytes, void *aux_data)
{
  (void)aux_data;
  void *dev_ptr;
  cudaError_t err = cudaMalloc(&dev_ptr, size_in_bytes);
  if (err) return NULL;
  //printf("Allocated device ptr %p\n", dev_ptr);
  return dev_ptr;
}

static axbStatus_t cuda_free(void *ptr_to_free, void *aux_data)
{
  (void)aux_data;
  cudaFree(ptr_to_free);
  return 0;
}


static axbStatus_t host_copyin(void *src, axbDataType_t src_type, void *dest, axbDataType_t dest_type, size_t n)
{
  if (src_type != AXB_REAL_DOUBLE || dest_type != AXB_REAL_DOUBLE) return 17590; // not yet supported

  //printf("Calling CUDA memcopy host to device %p\n", dest);
  return cudaMemcpy(dest, src, sizeof(double) * n, cudaMemcpyHostToDevice);
}

static axbStatus_t host_copyout(void *src, axbDataType_t src_type, void *dest, axbDataType_t dest_type, size_t n)
{
  if (src_type != AXB_REAL_DOUBLE || dest_type != AXB_REAL_DOUBLE) return 17590; // not yet supported
  //printf("Calling CUDA memcopy device %p to host\n", src);
  return cudaMemcpy(dest, src, sizeof(double) * n, cudaMemcpyDeviceToHost);
}

extern "C" axbStatus_t axbMemBackendRegister_CUDA(axbHandle_t handle)
{
  axbMemBackend_t cuda_backend;
  axbMemBackendCreate(&cuda_backend);

  // populate host_backend:
  axbMemBackendSetName(cuda_backend, "CUDA");
  axbMemBackendSetMalloc(cuda_backend, cuda_malloc);
  axbMemBackendSetFree(cuda_backend, cuda_free);

  axbMemBackendSetCopyIn(cuda_backend, host_copyin);
  axbMemBackendSetCopyOut(cuda_backend, host_copyout);

  // push into enclosing context identified by handle:
  axbMemBackendRegister(handle, cuda_backend);
  return 0;
}

