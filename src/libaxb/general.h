#ifndef LIBAXB_GENERAL_H_
#define LIBAXB_GENERAL_H_

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "libaxb/backend.h"


//#define AXB_ERRCHK(arg)     do { if (arg != 0) { fprintf(stderr, "ERROR in %s line %d: Code %d.\n", __FILE__, __LINE__, arg); assert( arg==0 ); return arg; } } while (0)
#define AXB_ERRCHK(arg)


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
