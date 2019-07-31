
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
  axbStatus_t status = axbOpBackendCreate(&cuda_backend); AXB_ERRCHK(status);

  // populate host_backend:
  status = axbOpBackendSetName(cuda_backend, "CUDA"); AXB_ERRCHK(status);

  axbOperationID_t op_id = 0;
  status = axbOpBackendAddOperation(cuda_backend, "vec-axpy", (axbStatus_t (*)(void))op_axpy, NULL, &op_id); AXB_ERRCHK(status);

  // push into enclosing context identified by handle:
  status = axbOpBackendRegister(handle, cuda_backend); AXB_ERRCHK(status);

  return 0;
}

