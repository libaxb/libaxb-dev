#ifndef LIBAXB_BACKEND_OPENCL_H_
#define LIBAXB_BACKEND_OPENCL_H_


#include "libaxb.h"

axbStatus_t axbMemBackendRegister_OpenCL(axbHandle_t handle);
axbStatus_t axbOpBackendRegister_OpenCL(axbHandle_t handle);

#ifdef LIBAXB_ENABLE_OPENCL

#include "CL/cl.h"

typedef struct axbMemOpenCL_s *axbMemOpenCL_t;

struct axbMemOpenCL_s
{
  cl_platform_id   platform;
  cl_context       context;
  cl_command_queue queue;

  cl_device_id     *devices;
  size_t           devices_size;
  size_t           devices_capacity;
};
#endif

#endif
