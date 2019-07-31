
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "libaxb.h"
#include "libaxb/backend.h"
#include "libaxb/general.h"

#include <cuda_runtime.h>

__global__
void kernel_axpy(int n, double *y, const double *alpha, const double *x)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) y[i] = *alpha * x[i] + y[i];
}


static axbStatus_t op_axpy(axbVec_t y, axbScalar_t alpha, axbVec_t x, void *aux_data)
{
  (void)aux_data;

  double *d_y     = (double*)y->data;
  double *d_alpha = (double*)alpha->data;
  double *d_x     = (double*)x->data;

  kernel_axpy<<<256, 256>>>((int)y->size, d_y, d_alpha, d_x);

  return 0;
}




extern "C" axbStatus_t axbOpBackendRegister_CUDA(axbHandle_t handle)
{
  axbOpBackend_t cuda_backend;
  axbOpBackendCreate(&cuda_backend);

  // populate host_backend:
  axbOpBackendSetName(cuda_backend, "CUDA");

  axbOperationID_t op_id = 0;
  axbOpBackendAddOperation(cuda_backend, "vec-axpy", (axbStatus_t (*)(void))op_axpy, NULL, &op_id);

  // push into enclosing context identified by handle:
  axbOpBackendRegister(handle, cuda_backend);

  return 0;
}

