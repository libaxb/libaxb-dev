
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "libaxb.h"
#include "libaxb/backend.h"
#include "libaxb/general.h"


axbStatus_t op_axpy(axbVec_t y, axbScalar_t alpha, axbVec_t x, void *aux_data)
{
  (void)aux_data;

  double *d_y     = y->data;
  double *d_alpha = (double*) alpha->data;
  double *d_x     = x->data;

  for (size_t i=0; i<y->size; ++i) d_y[i] += *d_alpha * d_x[i];

  return 0;
}




axbStatus_t axbOpBackendRegister_Host(axbHandle_t handle)
{
  axbOpBackend_t host_backend;
  axbOpBackendCreate(&host_backend);

  // populate host_backend:
  axbOpBackendSetName(host_backend, "host");

  axbOperationID_t op_id = 0;
  axbOpBackendAddOperation(host_backend, "vec-axpy", (axbStatus_t (*)(void))op_axpy, NULL, &op_id);

  // push into enclosing context identified by handle:
  axbOpBackendRegister(handle, host_backend);

  return 0;
}

