#ifndef LIBAXB_GENERAL_H_
#define LIBAXB_GENERAL_H_

#include <stdlib.h>
#include "libaxb/backend.h"

struct axbHandle_s
{
  int init;

  size_t            memBackends_size;
  size_t            memBackends_capacity;
  axbMemBackend_t  *memBackends;

  size_t            opBackends_size;
  size_t            opBackends_capacity;
  axbOpBackend_t   *opBackends;
};



struct axbScalar_s
{
  axbHandle_t handle;
  int init;

  size_t name_capacity;
  char *name;

  axbDataType_t datatype;
  void *data;

  axbMemBackend_t memBackend;
  axbOpBackend_t  opBackend;
};


struct axbVec_s
{
  axbHandle_t handle;
  int init;

  size_t size;

  size_t name_capacity;
  char *name;

  axbDataType_t datatype;
  void *data;

  axbMemBackend_t memBackend;
  axbOpBackend_t  opBackend;
};


#endif
