#include <stdlib.h>
#include <string.h>
#include "libaxb.h"
#include "libaxb/backend.h"

static void * host_malloc(size_t size_in_bytes, void *aux_data)
{
  (void)aux_data;
  return malloc(size_in_bytes);
}

static axbStatus_t host_free(void *ptr_to_free, void *aux_data)
{
  (void)aux_data;
  free(ptr_to_free);
  return 0;
}

axbStatus_t axbMemBackendRegister_Host(axbHandle_t handle)
{
  axbMemBackend_t host_backend;
  axbMemBackendCreate(&host_backend);

  // populate host_backend:
  axbMemBackendSetName(host_backend, "host");
  axbMemBackendSetMalloc(host_backend, host_malloc);
  axbMemBackendSetFree(host_backend, host_free);

  // push into enclosing context identified by handle:
  axbMemBackendRegister(handle, host_backend);
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

