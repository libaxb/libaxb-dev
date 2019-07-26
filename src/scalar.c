
#include <string.h>
#include "libaxb.h"
#include "libaxb/general.h"


axbStatus_t axbScalarCreateBegin(axbHandle_t handle, axbScalar_t *scalar)
{
  *scalar = malloc(sizeof(struct axbScalar_s));
  (*scalar)->handle = handle;
  (*scalar)->init = 896283;

  // set defaults:
  (*scalar)->name = malloc(10);
  (*scalar)->name[0] = 0;
  (*scalar)->name_capacity = 10;
  for (size_t i=0; i<sizeof((*scalar)->data); ++i) (*scalar)->data[0] = 0;

  (*scalar)->datatype = AXB_REAL_DOUBLE;
  (*scalar)->memBackend = handle->memBackends[0];
  (*scalar)->opBackend  = handle->opBackends[0];

  return 0;
}

axbStatus_t axbScalarSetDataType(axbScalar_t scalar, axbDataType_t datatype)
{
  scalar->datatype = datatype;
  return 0;
}
axbStatus_t axbScalarSetMemBackend(axbScalar_t scalar, axbMemBackend_t mem)
{
  scalar->memBackend = mem;
  return 0;
}
axbStatus_t axbScalarCreateEnd(axbScalar_t scalar)
{
  // nothing to do -> scalar already created!
  return 0;
}

axbStatus_t axbScalarSetValue(axbScalar_t scalar, void *value, axbDataType_t value_datatype)
{
  if (scalar->datatype != value_datatype) return 8713;

  double *val = (double *)value;
  double *d_data = (double *)scalar->data;
  *d_data = *val;
  return 0;
}

// convenience routine?
axbStatus_t axbScalarCreate(axbHandle_t handle, axbScalar_t *scalar, void *value, axbDataType_t datatype, axbMemBackend_t mem)
{
  axbScalarCreateBegin(handle, scalar);
  axbScalarSetDataType(*scalar, datatype);
  axbScalarSetMemBackend(*scalar, mem);
  axbScalarCreateEnd(*scalar);
  axbScalarSetValue(*scalar, value, datatype);
  return 0;
}


axbStatus_t axbScalarDestroy(axbScalar_t scalar)
{
  if (scalar->init != 896283) return 896283;
  scalar->init += 1;

  free(scalar->name);
  free(scalar);
  return 0;
}
