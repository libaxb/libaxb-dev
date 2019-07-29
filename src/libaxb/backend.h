#ifndef LIBAXB_BACKEND_H_
#define LIBAXB_BACKEND_H_

#include "libaxb.h"

struct axbMemBackend_s
{
  size_t name_capacity;
  char *name;

  void *      (*op_malloc)(size_t, void*);
  axbStatus_t (*op_free)  (void*, void*);

  void *impl;   //pimpl idiom for backends to drop in their specific stuff
};


////////////////////

struct axbOpDescriptor_s
{
  axbStatus_t (*func)(void);   // type-agnostic function pointer
  void  *func_data;

  char  name[32];

  axbOperationID_t id;
};
typedef struct axbOpDescriptor_s axbOpDescriptor_t;

struct axbOpBackend_s
{
  size_t name_capacity;
  char *name;

  // table of operations:
  size_t op_table_size;
  size_t op_table_capacity;
  axbOpDescriptor_t *op_table;

  // more operations to be added
  void *impl;
};


axbStatus_t axbMemBackendRegisterDefaults(axbHandle_t handle);
axbStatus_t axbOpBackendRegisterDefaults(axbHandle_t handle);

#endif
