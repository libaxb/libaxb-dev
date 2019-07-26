
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

axbStatus_t axbOpBackendSetName(axbOpBackend_t ops, const char *name)
{
  size_t len = strlen(name);
  if (len > ops->name_capacity) {
    free(ops->name);
    ops->name_capacity = len + 2;
    ops->name = malloc(ops->name_capacity);
  }
  for (size_t i=0; i<=len; ++i) ops->name[i] = name[i];

  return 0;
}

axbStatus_t axbOpBackendGetName(axbOpBackend_t ops, const char **name)
{
  *name = ops->name;
  return 0;
}


axbStatus_t axbOpBackendRegister_Host(axbHandle_t handle)
{
  axbOpBackend_t host_backend;
  axbOpBackendCreate(&host_backend);

  // populate host_backend:
  axbOpBackendSetName(host_backend, "host");

  axbOpBackendSetOperation(host_backend, 0, (void (*)(void)) op_axpy, NULL);

  // push into enclosing context identified by handle:
  axbOpBackendRegister(handle, host_backend);

  return 0;
}


axbStatus_t axbOpBackendSetOperation(axbOpBackend_t ops, const axbOperationID_t op_id, void (*func)(void), void *aux_data)
{
  // TODO: Dispatch with respect to `op_id`
  (void)op_id;
  ops->op_axpy = (axbStatus_t (*)(axbVec_t, axbScalar_t, axbVec_t, void*))func;
  ops->op_axpy_data = aux_data;
  return 0;
}
