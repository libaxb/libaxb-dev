
#include <stdlib.h>
#include "libaxb.h"

#include "libaxb/general.h"


axbStatus_t axbMemBackendRegister(axbHandle_t handle, axbMemBackend_t mem)
{
  // check if there is still room for another backend. If not, resize array:
  if (handle->memBackends_size == handle->memBackends_capacity) {
    axbMemBackend_t *old_array = handle->memBackends;
    handle->memBackends_capacity *= 2;
    handle->memBackends = malloc(handle->memBackends_capacity * sizeof(axbMemBackend_t));

    // copy over existing backends:
    for (size_t i=0; i<handle->memBackends_size; ++i) handle->memBackends[i] = old_array[i];
    free(old_array);
  }

  handle->memBackends[handle->memBackends_size] = mem;
  handle->memBackends_size += 1;
  return 0;
}

axbStatus_t axbMemBackendGetAll(axbHandle_t handle, axbMemBackend_t **mem, int *mem_size)
{
  *mem = handle->memBackends;
  *mem_size = handle->memBackends_size;
  return 0;
}



///////////////

axbStatus_t axbOpBackendRegister(axbHandle_t handle, axbOpBackend_t mem)
{
  // check if there is still room for another backend. If not, resize array:
  if (handle->opBackends_size == handle->opBackends_capacity) {
    axbOpBackend_t *old_array = handle->opBackends;
    handle->opBackends_capacity *= 2;
    handle->opBackends = malloc(handle->opBackends_capacity * sizeof(axbOpBackend_t));

    // copy over existing backends:
    for (size_t i=0; i<handle->opBackends_size; ++i) handle->opBackends[i] = old_array[i];
    free(old_array);
  }

  handle->opBackends[handle->opBackends_size] = mem;
  handle->opBackends_size += 1;
  return 0;
}

axbStatus_t axbOpBackendGetAll(axbHandle_t handle, axbOpBackend_t **ops, int *ops_size)
{
  *ops = handle->opBackends;
  *ops_size = handle->opBackends_size;
  return 0;
}


/////////////////

axbStatus_t axbInit(axbHandle_t *handle)
{
  *handle = malloc(sizeof(struct axbHandle_s));
  (*handle)->init = 424242; // magic number to indicate proper initialization (and track repeated free's)

  (*handle)->memBackends_capacity = 10;
  (*handle)->memBackends_size = 0;
  (*handle)->memBackends = malloc((*handle)->memBackends_capacity * sizeof(axbMemBackend_t));

  axbMemBackendRegisterDefaults(*handle);

  (*handle)->opBackends_capacity = 10;
  (*handle)->opBackends_size = 0;
  (*handle)->opBackends = malloc((*handle)->opBackends_capacity * sizeof(axbMemBackend_t));

  axbOpBackendRegisterDefaults(*handle);

  return 0;
}



axbStatus_t axbFinalize(axbHandle_t handle)
{
  // guard against multiple free()
  if (handle->init != 424242) return 424242;
  handle->init += 1; // helper to track multiple free's on handle

  // free mem backends
  for (size_t i=0; i<handle->memBackends_size; ++i) {
    axbMemBackendDestroy(handle->memBackends[i]);
  }
  free(handle->memBackends);

  // free op backends
  for (size_t i=0; i<handle->opBackends_size; ++i) {
    axbOpBackendDestroy(handle->opBackends[i]);
  }
  free(handle->opBackends);

  free(handle);
  return 0;
}

