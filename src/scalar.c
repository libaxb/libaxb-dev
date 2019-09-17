
#include <string.h>
#include "libaxb.h"
#include "libaxb/general.h"


axbStatus_t axbScalarCreateBegin(struct axbHandle_s *handle, struct axbScalar_s **scalar)
{
  *scalar = malloc(sizeof(struct axbScalar_s));
  (*scalar)->handle = handle;
  (*scalar)->init = 896283;

  // set defaults:
  (*scalar)->name = malloc(10);
  (*scalar)->name[0] = 0;
  (*scalar)->name_capacity = 10;

  (*scalar)->datatype = AXB_REAL_DOUBLE;
  (*scalar)->memBackend = handle->memBackends[0];
  (*scalar)->opBackend  = handle->opBackends[0];

  return 0;
}

axbStatus_t axbScalarSetDataType(struct axbScalar_s *scalar, axbDataType_t datatype)
{
  scalar->datatype = datatype;
  return 0;
}

axbStatus_t axbScalarSetMemBackend(struct axbScalar_s *scalar, struct axbMemBackend_s *mem)
{
  if (mem) scalar->memBackend = mem;
  return 0;
}

axbStatus_t axbScalarCreateEnd(struct axbScalar_s *scalar)
{
  return scalar->memBackend->op_malloc( &(scalar->data), sizeof(double), scalar->memBackend->impl);
}

axbStatus_t axbScalarSetValue(struct axbScalar_s *scalar, void *value, axbDataType_t value_datatype)
{
  return axbMemBackendCopyIn(scalar->memBackend, value, value_datatype, scalar->data, scalar->datatype, 1);
}

axbStatus_t axbScalarGetValue(const struct axbScalar_s *scalar, void *value, axbDataType_t value_datatype)
{
  return axbMemBackendCopyOut(scalar->memBackend, scalar->data, scalar->datatype, value, value_datatype, 1);
}

// convenience routine?
axbStatus_t axbScalarCreate(struct axbHandle_s *handle, struct axbScalar_s **scalar, void *value, axbDataType_t datatype, struct axbMemBackend_s *mem)
{
  axbScalarCreateBegin(handle, scalar);
  axbScalarSetDataType(*scalar, datatype);
  axbScalarSetMemBackend(*scalar, mem);
  axbScalarCreateEnd(*scalar);
  axbScalarSetValue(*scalar, value, datatype);
  return 0;
}


axbStatus_t axbScalarDestroy(struct axbScalar_s *scalar)
{
  if (scalar->init != 896283) return 896283;
  scalar->init += 1;

  scalar->memBackend->op_free(scalar->data, NULL);
  free(scalar->name);
  free(scalar);
  return 0;
}
