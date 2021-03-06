#ifndef LIBAXB_BACKEND_H_
#define LIBAXB_BACKEND_H_

#include "libaxb.h"


struct axbMemBackend_s
{
  size_t name_capacity;
  char *name;

  axbStatus_t (*op_malloc)(void**, size_t, void*);
  axbStatus_t (*op_free)  (void*, void*);

  axbStatus_t (*op_copyin) (void*, axbDataType_t, void*, axbDataType_t, size_t, void *);
  axbStatus_t (*op_copyout)(void*, axbDataType_t, void*, axbDataType_t, size_t, void *);

  axbStatus_t (*destroy)(void*);

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

  axbStatus_t (*destroy)(void*);

  // more operations to be added
  void *impl;
};



#ifdef __cplusplus
extern "C"
{
#endif

axbStatus_t axbMemBackendRegisterDefaults(struct axbHandle_s *handle);
axbStatus_t axbMemBackendSetValues(struct axbMemBackend_s *mem, void *src, axbDataType_t src_data, void *dest, void *dest_data);
axbStatus_t axbMemBackendGetValues(struct axbMemBackend_s *mem, void *src, axbDataType_t src_data, void *dest, void *dest_data);

axbStatus_t axbOpBackendRegisterDefaults(struct axbHandle_s *handle);


#ifdef __cplusplus
} // extern "C"
#endif


#endif
