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
