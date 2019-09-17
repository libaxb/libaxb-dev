
#include <stdlib.h>
#include <string.h>
#include "libaxb.h"

#include "libaxb/general.h"


axbStatus_t axbMemBackendRegister(struct axbHandle_s *handle, struct axbMemBackend_s *mem)
{
  // check if there is still room for another backend. If not, resize array:
  if (handle->memBackends_size == handle->memBackends_capacity) {
    struct axbMemBackend_s **old_array = handle->memBackends;
    handle->memBackends_capacity *= 2;
    handle->memBackends = malloc(handle->memBackends_capacity * sizeof(struct axbMemBackend_s *));

    // copy over existing backends:
    for (size_t i=0; i<handle->memBackends_size; ++i) handle->memBackends[i] = old_array[i];
    free(old_array);
  }

  handle->memBackends[handle->memBackends_size] = mem;
  handle->memBackends_size += 1;
  return 0;
}

axbStatus_t axbMemBackendGetAll(struct axbHandle_s *handle, struct axbMemBackend_s ***mem, size_t *mem_size)
{
  *mem = handle->memBackends;
  *mem_size = handle->memBackends_size;
  return 0;
}

axbStatus_t axbMemBackendGetByName(struct axbHandle_s *handle, struct axbMemBackend_s **mem, const char *name)
{
  *mem = NULL;
  for (size_t i=0; i<handle->memBackends_size; ++i){
    if ( strcmp(handle->memBackends[i]->name, name) == 0) {
      *mem = handle->memBackends[i];
      return 0;
    }
  }
  return 0;
}


///////////////

axbStatus_t axbOpBackendRegister(struct axbHandle_s *handle, struct axbOpBackend_s *ops)
{
  // check if there is still room for another backend. If not, resize array:
  if (handle->opBackends_size == handle->opBackends_capacity) {
    struct axbOpBackend_s **old_array = handle->opBackends;
    handle->opBackends_capacity *= 2;
    handle->opBackends = malloc(handle->opBackends_capacity * sizeof(struct axbOpBackend_s *));

    // copy over existing backends:
    for (size_t i=0; i<handle->opBackends_size; ++i) handle->opBackends[i] = old_array[i];
    free(old_array);
  }

  handle->opBackends[handle->opBackends_size] = ops;
  handle->opBackends_size += 1;
  return 0;
}

axbStatus_t axbOpBackendGetAll(struct axbHandle_s *handle, struct axbOpBackend_s ***ops, size_t *ops_size)
{
  *ops = handle->opBackends;
  *ops_size = handle->opBackends_size;
  return 0;
}

axbStatus_t axbOpBackendGetByName(struct axbHandle_s *handle, struct axbOpBackend_s **ops, const char *name)
{
  *ops = NULL;
  for (size_t i=0; i<handle->opBackends_size; ++i){
    if ( strcmp(handle->opBackends[i]->name, name) == 0) {
      *ops = handle->opBackends[i];
      return 0;
    }
  }
  return 0;
}

/////////////////

axbStatus_t axbInit(struct axbHandle_s **handle)
{
  *handle = malloc(sizeof(struct axbHandle_s));
  (*handle)->init = 424242; // magic number to indicate proper initialization (and track repeated free's)

  (*handle)->memBackends_capacity = 10;
  (*handle)->memBackends_size = 0;
  (*handle)->memBackends = malloc((*handle)->memBackends_capacity * sizeof(struct axbMemBackend_s *));

  axbMemBackendRegisterDefaults(*handle);

  (*handle)->opBackends_capacity = 10;
  (*handle)->opBackends_size = 0;
  (*handle)->opBackends = malloc((*handle)->opBackends_capacity * sizeof(struct axbMemBackend_s *));

  axbOpBackendRegisterDefaults(*handle);

  return 0;
}



axbStatus_t axbFinalize(struct axbHandle_s *handle)
{
  // guard against multiple free()
  if (handle->init != 424242) return 424242;
  handle->init += 1; // helper to track multiple free's on handle

  // free op backends
  // Note: op-backends may depend on mem-backends, so it's important to free the op-backends first.
  for (size_t i=0; i<handle->opBackends_size; ++i) {
    axbOpBackendDestroy(handle->opBackends[i]);
  }
  free(handle->opBackends);

  // free mem backends
  for (size_t i=0; i<handle->memBackends_size; ++i) {
    axbMemBackendDestroy(handle->memBackends[i]);
  }
  free(handle->memBackends);

  free(handle);
  return 0;
}

