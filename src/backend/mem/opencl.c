#ifdef LIBAXB_ENABLE_OPENCL



#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "libaxb.h"
#include "libaxb/general.h"
#include "libaxb/backend.h"
#include "libaxb/backend/opencl.h"


static axbStatus_t opencl_malloc(void **ptr, size_t size_in_bytes, void *aux_data)
{
  axbMemOpenCL_t opencl_context = aux_data;

  cl_int status;
  cl_mem mem = clCreateBuffer(opencl_context->context, CL_MEM_READ_WRITE, size_in_bytes, NULL, &status); AXB_ERRCHK(status);
  *ptr = mem;
  return status;
}

static axbStatus_t opencl_free(void *ptr_to_free, void *aux_data)
{
  (void)aux_data;
  cl_mem mem = ptr_to_free;
  cl_int status = clReleaseMemObject(mem); AXB_ERRCHK(status);
  return 0;
}

static axbStatus_t opencl_copyin(void *src, axbDataType_t src_type, void *dest, axbDataType_t dest_type, size_t n, void *aux_data)
{
  if (src_type != AXB_REAL_DOUBLE || dest_type != AXB_REAL_DOUBLE) return 16590; // not yet supported

  axbMemOpenCL_t opencl_context = aux_data;
  cl_mem cl_dest = dest;
  cl_int status = clEnqueueWriteBuffer(opencl_context->queue, cl_dest, CL_TRUE, 0, sizeof(double) * n, src, 0, NULL, NULL); AXB_ERRCHK(status);

  return status;
}

static axbStatus_t opencl_copyout(void *src, axbDataType_t src_type, void *dest, axbDataType_t dest_type, size_t n, void *aux_data)
{
  if (src_type != AXB_REAL_DOUBLE || dest_type != AXB_REAL_DOUBLE) return 16590; // not yet supported

  axbMemOpenCL_t opencl_context = aux_data;
  cl_mem cl_src = src;
  cl_int status = clEnqueueReadBuffer(opencl_context->queue, cl_src, CL_TRUE, 0, sizeof(double) * n, dest, 0, NULL, NULL); AXB_ERRCHK(status);

  return status;
}


static axbStatus_t destroyOpenCLContext(void *impl)
{
  axbMemOpenCL_t opencl_context = (axbMemOpenCL_t)impl;
  cl_int status = clReleaseCommandQueue(opencl_context->queue);  AXB_ERRCHK(status);
  status = clReleaseContext(opencl_context->context);  AXB_ERRCHK(status);
  free(opencl_context->devices);
  free(opencl_context);
  return 0;
}


axbStatus_t axbMemBackendCreate_OpenCL(axbMemBackend_t opencl_mem_backend)
{
  axbMemOpenCL_t opencl_context = malloc(sizeof(struct axbMemOpenCL_s));

  // set function pointers:
  opencl_mem_backend->impl = opencl_context;

  // set up platform:
  cl_uint platforms_capacity = 42;
  cl_uint platforms_size = 0;
  cl_platform_id *platforms = malloc(platforms_capacity*sizeof(cl_platform_id));
  cl_int status = clGetPlatformIDs(platforms_capacity, platforms, &platforms_size); AXB_ERRCHK(status);
  opencl_context->platform = platforms[0];

  // get devices for context:
  cl_device_type device_type = CL_DEVICE_TYPE_ALL;
  opencl_context->devices_capacity = 42;
  opencl_context->devices_size = 0;
  opencl_context->devices = malloc(opencl_context->devices_capacity*sizeof(cl_device_id));
  cl_uint num_devices, capacity = (cl_uint)opencl_context->devices_capacity; // OpenCL expects cl_uints...
  status = clGetDeviceIDs(opencl_context->platform, device_type, capacity, opencl_context->devices, &num_devices); AXB_ERRCHK(status);
  opencl_context->devices_size = (size_t)num_devices;

  // set up context:
  cl_context_properties context_props[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)opencl_context->platform, 0 };
  opencl_context->context = clCreateContext(context_props, num_devices, opencl_context->devices, NULL, NULL, &status); AXB_ERRCHK(status);

  // add a command queue:
  opencl_context->queue = clCreateCommandQueue(opencl_context->context, opencl_context->devices[0], 0, &status); AXB_ERRCHK(status);

  // clean up temporary arrays:
  free(platforms);


  return 0;
}


axbStatus_t axbMemBackendRegister_OpenCL(axbHandle_t handle)
{
  axbMemBackend_t opencl_backend;
  axbStatus_t status = axbMemBackendCreate(&opencl_backend); AXB_ERRCHK(status);

  // spawn OpenCL context:
  status = axbMemBackendCreate_OpenCL(opencl_backend); AXB_ERRCHK(status);

  // populate function pointers:
  status = axbMemBackendSetName(opencl_backend, "OpenCL"); AXB_ERRCHK(status);
  status = axbMemBackendSetMalloc(opencl_backend, opencl_malloc); AXB_ERRCHK(status);
  status = axbMemBackendSetFree(opencl_backend, opencl_free); AXB_ERRCHK(status);

  status = axbMemBackendSetCopyIn(opencl_backend, opencl_copyin); AXB_ERRCHK(status);
  status = axbMemBackendSetCopyOut(opencl_backend, opencl_copyout); AXB_ERRCHK(status);

  status = axbMemBackendSetDestroy(opencl_backend, destroyOpenCLContext); AXB_ERRCHK(status);

  // push into enclosing context identified by handle:
  status = axbMemBackendRegister(handle, opencl_backend); AXB_ERRCHK(status);
  return 0;
}

#else

#include "libaxb.h"

axbStatus_t axbMemBackendRegister_OpenCL(axbHandle_t handle)
{
  return 0;
}

#endif
