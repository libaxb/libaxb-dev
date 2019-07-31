
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "libaxb.h"
#include "libaxb/backend.h"
#include "libaxb/general.h"


static axbStatus_t op_axpy(axbVec_t y, axbScalar_t alpha, axbVec_t x, void *aux_data)
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
  axbStatus_t status = axbOpBackendCreate(&host_backend); AXB_ERRCHK(status);

  // populate host_backend:
  status = axbOpBackendSetName(host_backend, "host"); AXB_ERRCHK(status);

  axbOperationID_t op_id = 0;
  status = axbOpBackendAddOperation(host_backend, "vec-axpy", (axbStatus_t (*)(void))op_axpy, NULL, &op_id); AXB_ERRCHK(status);

  // push into enclosing context identified by handle:
  status = axbOpBackendRegister(handle, host_backend); AXB_ERRCHK(status);

  return 0;
}

