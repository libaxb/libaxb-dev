
#include <string.h>
#include "libaxb.h"
#include "libaxb/general.h"

axbStatus_t axbVecCreateBegin(axbHandle_t handle, axbVec_t *vec)
{
  *vec = malloc(sizeof(struct axbVec_s));
  (*vec)->handle = handle;
  (*vec)->init = 119346;

  // set defaults:
  (*vec)->name = malloc(10);
  (*vec)->name[0] = 0;
  (*vec)->name_capacity = 10;
  (*vec)->data = NULL;

  (*vec)->datatype = AXB_REAL_DOUBLE;
  (*vec)->memBackend = handle->memBackends[0];
  (*vec)->opBackend  = handle->opBackends[0];

  return 0;
}

axbStatus_t axbVecSetSize(axbVec_t vec, size_t size)
{
  // TODO: Check that VecCreateEnd() has not been called yet!
  vec->size = size;
  return 0;
}
axbStatus_t axbVecGetSize(axbVec_t vec, size_t *size)
{
  *size = vec->size;
  return 0;
}

axbStatus_t axbVecSetDataType(axbVec_t vec, axbDataType_t datatype)
{
  // TODO: Check that VecCreateEnd() has not been called yet!
  vec->datatype = datatype;
  return 0;
}
axbStatus_t axbVecGetDataType(axbVec_t vec, axbDataType_t *datatype)
{
  *datatype = vec->datatype;
  return 0;
}

axbStatus_t axbVecSetMemBackend(axbVec_t vec, axbMemBackend_t backend)
{
  vec->memBackend = backend;
  return 0;
}
axbStatus_t axbVecGetMemBackend(axbVec_t vec, axbMemBackend_t *backend)
{
  *backend = vec->memBackend;
  return 0;
}

axbStatus_t axbVecSetOpBackend(axbVec_t vec, axbOpBackend_t backend)
{
  vec->opBackend = backend;
  return 0;
}
axbStatus_t axbVecGetOpBackend(axbVec_t vec, axbOpBackend_t *backend)
{
  *backend = vec->opBackend;
  return 0;
}


axbStatus_t axbVecCreateEnd(axbVec_t vec)
{
  if (vec->size > 0) {
    vec->data = vec->memBackend->op_malloc(sizeof(double) * vec->size, vec->memBackend->impl);
    if (!vec->data) return 123; // out of memory
  }
  return 0;
}

axbStatus_t axbVecSetName(axbVec_t vec, const char *name)
{
  size_t len = strlen(name);
  if (len > vec->name_capacity) {
    free(vec->name);
    vec->name = malloc(len + 2);
  }
  for (size_t i=0; i<=len; ++i) vec->name[i] = name[i];
  return 0;
}
axbStatus_t axbVecGetName(axbVec_t vec, const char **name)
{
  *name = vec->name;
  return 0;
}

axbStatus_t axbVecSetValues(axbVec_t vec, void *values, axbDataType_t values_datatype)
{
  if (values_datatype != vec->datatype) return 1237;  // Data conversion not yet implemented!

  double *d_data = vec->data;
  double *d_values = values;
  for (size_t i=0; i<vec->size; ++i)
  {
    d_data[i] = d_values[i];
  }
  return 0;
}
axbStatus_t axbVecGetValues(axbVec_t vec, void *values, axbDataType_t values_datatype)
{
  if (values_datatype != vec->datatype) return 1237;  // Data conversion not yet implemented!

  double *d_data = vec->data;
  double *d_values = values;
  for (size_t i=0; i<vec->size; ++i)
  {
    d_values[i] = d_data[i];
  }

  return 0;
}


// operations

/** @brief y = alpha * x + y */
axbStatus_t axbVecAXPY(axbVec_t y, axbScalar_t alpha, axbVec_t x)
{
  y->opBackend->op_axpy(y, alpha, x, y->opBackend->impl);
}
// more to follow

axbStatus_t axbVecDestroy(axbVec_t vec)
{
  if (vec->init != 119346) return 119346;
  vec->init += 1;

  free(vec->data);
  free(vec->name);
  free(vec);
  return 0;
}
