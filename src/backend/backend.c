
#include <string.h>

#include "libaxb.h"
#include "libaxb/general.h"
#include "libaxb/backend/host.h"
#include "libaxb/backend/cuda.h"
#include "libaxb/backend/opencl.h"
#include "libaxb/backend.h"

axbStatus_t axbMemBackendCreate(axbMemBackend_t *mem)
{
  *mem = malloc(sizeof(struct axbMemBackend_s));

  (*mem)->name_capacity = 10;
  (*mem)->name = malloc((*mem)->name_capacity);
  (*mem)->impl = NULL;
  (*mem)->op_malloc  = NULL;
  (*mem)->op_free    = NULL;
  (*mem)->op_copyin  = NULL;
  (*mem)->op_copyout = NULL;
  (*mem)->destroy    = NULL;

  return 0;
}

axbStatus_t axbMemBackendRegisterDefaults(axbHandle_t handle)
{
  axbMemBackendRegister_Host(handle);
  axbMemBackendRegister_CUDA(handle);
  axbMemBackendRegister_OpenCL(handle);
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

axbStatus_t axbMemBackendSetMalloc(axbMemBackend_t mem, axbStatus_t (*func)(void **, size_t, void *))
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
  return mem->op_malloc(ptr, num_bytes, mem->impl);
}

axbStatus_t axbMemBackendFree(axbMemBackend_t mem, void *ptr)
{
  return mem->op_free(ptr, mem->impl);
}


axbStatus_t axbMemBackendSetCopyIn(axbMemBackend_t mem, axbStatus_t (*func)(void *, axbDataType_t, void *, axbDataType_t, size_t, void *))
{
  mem->op_copyin = func;
  return 0;
}
axbStatus_t axbMemBackendSetCopyOut(axbMemBackend_t mem, axbStatus_t (*func)(void *, axbDataType_t, void *, axbDataType_t, size_t, void *))
{
  mem->op_copyout = func;
  return 0;
}

axbStatus_t axbMemBackendCopyIn(axbMemBackend_t mem, void *src, axbDataType_t src_type, void *dest, axbDataType_t dest_type, size_t n)
{
  return mem->op_copyin(src, src_type, dest, dest_type, n, mem->impl);
}
axbStatus_t axbMemBackendCopyOut(axbMemBackend_t mem, void *src, axbDataType_t src_type, void *dest, axbDataType_t dest_type, size_t n)
{
  return mem->op_copyout(src, src_type, dest, dest_type, n, mem->impl);
}

axbStatus_t axbMemBackendSetDestroy(axbMemBackend_t mem, axbStatus_t (*func)(void*))
{
  mem->destroy = func;
  return 0;
}

axbStatus_t axbMemBackendDestroy(axbMemBackend_t mem)
{
  if (mem->destroy) {
    axbStatus_t status = mem->destroy(mem->impl); AXB_ERRCHK(status);
  }
  free(mem->name);
  free(mem);
  return 0;
}


////////////////


axbStatus_t axbOpBackendCreate(axbOpBackend_t *ops)
{
  *ops = malloc(sizeof(struct axbOpBackend_s));

  (*ops)->name = NULL;
  (*ops)->name_capacity = 0;
  (*ops)->impl = NULL;
  (*ops)->op_table_capacity = 255;
  (*ops)->op_table_size = 0;
  (*ops)->op_table = malloc((*ops)->op_table_capacity*sizeof(axbOpDescriptor_t));

  (*ops)->destroy = NULL;

  return 0;
}

axbStatus_t axbOpBackendSetName(axbOpBackend_t ops, const char *name)
{
  size_t len = strlen(name);
  if (!ops->name)
  {
    ops->name_capacity = len + 2;
    ops->name = malloc(ops->name_capacity);
  }

  if (len > sizeof(ops->name)) return 9753;  // op name too large
  for (size_t i=0; i<=len; ++i) ops->name[i] = name[i];

  return 0;
}

axbStatus_t axbOpBackendGetName(axbOpBackend_t ops, const char **name)
{
  *name = ops->name;
  return 0;
}

axbStatus_t axbOpBackendSetDestroy(axbOpBackend_t ops, axbStatus_t (*func)(void*))
{
  ops->destroy = func;
  return 0;
}

axbStatus_t axbOpBackendDestroy(axbOpBackend_t ops)
{
  if (ops->destroy) {
    axbStatus_t status = ops->destroy(ops->impl); AXB_ERRCHK(status);
  }
  free(ops->op_table);
  free(ops->name);
  free(ops);
  return 0;
}


axbStatus_t axbOpBackendAddOperation(axbOpBackend_t ops, const char *op_name, axbStatus_t (*op_func)(void), void *op_data, axbOperationID_t *op_id)
{
  if (ops->op_table_size == ops->op_table_capacity) {  // Expand op-table if already fully filled up
    axbOpDescriptor_t *old_op_table = ops->op_table;
    ops->op_table_capacity *= 2;
    ops->op_table = malloc(ops->op_table_capacity * sizeof(axbOpDescriptor_t));
    for (size_t i=0; i<ops->op_table_size; ++i) ops->op_table[i] = old_op_table[i];
  }

  axbOpDescriptor_t *op = ops->op_table + ops->op_table_size;
  op->func = op_func;
  op->func_data = op_data;
  op->name[sizeof(op->name)-1] = 0; // make sure string is terminated
  strncpy(op->name, op_name, sizeof(op->name)-1);
  op->id = ops->op_table_size;
  *op_id = op->id;
  ops->op_table_size += 1;
  return 0;
}


axbStatus_t axbOpBackendRegisterDefaults(axbHandle_t handle)
{
  axbStatus_t status;
  status = axbOpBackendRegister_Host(handle); AXB_ERRCHK(status);
  status = axbOpBackendRegister_CUDA(handle); AXB_ERRCHK(status);
  status = axbOpBackendRegister_OpenCL(handle); AXB_ERRCHK(status);
  return 0;
}


////////////////////


