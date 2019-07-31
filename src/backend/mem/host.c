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

static axbStatus_t host_copyin(void *src, axbDataType_t src_type, void *dest, axbDataType_t dest_type, size_t n, void *aux_data)
{
  if (src_type != AXB_REAL_DOUBLE || dest_type != AXB_REAL_DOUBLE) return 16590; // not yet supported

  double *d_src  = src;
  double *d_dest = dest;
  for (size_t i=0; i<n; ++i)
  {
    d_dest[i] = d_src[i];
  }

  (void)aux_data;
  return 0;
}

static axbStatus_t host_copyout(void *src, axbDataType_t src_type, void *dest, axbDataType_t dest_type, size_t n, void *aux_data)
{
  if (src_type != AXB_REAL_DOUBLE || dest_type != AXB_REAL_DOUBLE) return 16590; // not yet supported

  double *d_src  = src;
  double *d_dest = dest;
  for (size_t i=0; i<n; ++i)
  {
    d_dest[i] = d_src[i];
  }

  (void)aux_data;
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

  axbMemBackendSetCopyIn(host_backend, host_copyin);
  axbMemBackendSetCopyOut(host_backend, host_copyout);

  // push into enclosing context identified by handle:
  axbMemBackendRegister(handle, host_backend);
  return 0;
}
