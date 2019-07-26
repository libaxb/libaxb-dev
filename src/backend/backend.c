

#include "libaxb.h"

#include "libaxb/backend/host.h"
#include "libaxb/backend.h"

axbStatus_t axbMemBackendCreate(axbMemBackend_t *mem)
{
  *mem = malloc(sizeof(struct axbMemBackend_s));

  (*mem)->name_capacity = 10;
  (*mem)->name = malloc((*mem)->name_capacity);
  (*mem)->impl = NULL;
  (*mem)->op_malloc = NULL;
  (*mem)->op_free = NULL;

  return 0;
}

axbStatus_t axbMemBackendDestroy(axbMemBackend_t mem)
{
  free(mem->name);
  free(mem);
  return 0;
}

axbStatus_t axbMemBackendRegisterDefaults(axbHandle_t handle)
{
  axbMemBackendRegister_Host(handle);
  return 0;
}



////////////////


axbStatus_t axbOpBackendCreate(axbOpBackend_t *ops)
{
  *ops = malloc(sizeof(struct axbOpBackend_s));

  (*ops)->name_capacity = 10;
  (*ops)->name = malloc((*ops)->name_capacity);
  (*ops)->impl = NULL;
  (*ops)->op_axpy = NULL;
  (*ops)->op_axpy_data = NULL;

  return 0;
}

axbStatus_t axbOpBackendDestroy(axbOpBackend_t ops)
{
  free(ops->name);
  free(ops);
  return 0;
}

axbStatus_t axbOpBackendRegisterDefaults(axbHandle_t handle)
{
  axbOpBackendRegister_Host(handle);
  return 0;
}


////////////////////


