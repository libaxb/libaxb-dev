
#include <string.h>

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


axbStatus_t axbMemBackendSetName(axbMemBackend_t ops, const char *name)
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

axbStatus_t axbMemBackendGetName(axbMemBackend_t ops, const char **name)
{
  *name = ops->name;
  return 0;
}

axbStatus_t axbMemBackendSetMalloc(axbMemBackend_t mem, void *(*func)(size_t, void *))
{
  mem->op_malloc = func;
  return 0;
}
axbStatus_t axbMemBackendSetFree(axbMemBackend_t mem, axbStatus_t (*func)(void *, void *))
{
  mem->op_free = func;
  return 0;
}


axbStatus_t axbMemBackendMalloc(axbMemBackend_t mem, size_t num_bytes, void **ptr)
{
  *ptr = mem->op_malloc(num_bytes, mem->impl);
  return 0;
}

axbStatus_t axbMemBackendFree(axbMemBackend_t mem, void *ptr)
{
  mem->op_free(ptr, mem->impl);
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

axbStatus_t axbOpBackendDestroy(axbOpBackend_t ops)
{
  free(ops->name);
  free(ops);
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


axbStatus_t axbOpBackendRegisterDefaults(axbHandle_t handle)
{
  axbOpBackendRegister_Host(handle);
  return 0;
}


////////////////////


