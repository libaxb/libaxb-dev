#ifndef LIBAXB_BACKEND_H_
#define LIBAXB_BACKEND_H_

#include "libaxb.h"

#include "libaxb/backend.h"

struct axbMemBackend_s
{
  size_t name_capacity;
  char *name;

  void *      (*op_malloc)(size_t, void*);
  axbStatus_t (*op_free)  (void*, void*);

  void *impl;   //pimpl idiom for backends to drop in their specific stuff
};


struct axbOpBackend_s
{
  size_t name_capacity;
  char *name;

  // table of operations:
  axbStatus_t (*op_axpy)(axbVec_t y, axbScalar_t alpha, axbVec_t x, void *data);
  void *op_axpy_data;

  // more operations to be added
  void *impl;
};


axbStatus_t axbMemBackendRegisterDefaults(axbHandle_t handle);
axbStatus_t axbOpBackendRegisterDefaults(axbHandle_t handle);

#endif
