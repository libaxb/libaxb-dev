
#ifdef LIBAXB_ENABLE_OPENCL

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "libaxb.h"
#include "libaxb/backend.h"
#include "libaxb/general.h"
#include "libaxb/backend/opencl.h"

typedef struct axbOpOpenCL_s *axbOpOpenCL_t;

struct axbOpOpenCL_s
{
  // from struct axbMemOpenCL_s:
  // TODO: Think about referencing axbMemOpenCL_s instead
  cl_platform_id   platform;
  cl_context       context;
  cl_command_queue queue;

  // additional members for operations:
  cl_program *programs;
  size_t     programs_size;     // number of actual programs
  size_t     programs_capacity; // allocated size of programs-buffer

  cl_kernel  *kernels;
  size_t     kernels_size;      // number of actual kernels
  size_t     kernels_capacity;  // allocated size of kernels-buffer
};



static const char * my_compute_program =
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable \n"
"__kernel void axpy(\n"
"          __global       double *y,\n"
"          __global const double *alpha, \n"
"          __global const double *x,\n"
"          unsigned int size) \n"
"{ \n"
"  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0))\n"
"    y[i] += *alpha * x[i];\n"
"};\n\n";



static axbStatus_t op_axpy_opencl(axbVec_t y, axbScalar_t alpha, axbVec_t x, void *aux_data)
{
  cl_int status;
  axbOpOpenCL_t op_opencl = (axbOpOpenCL_t)aux_data;

  cl_mem cl_y     = y->data;
  cl_mem cl_alpha = alpha->data;
  cl_mem cl_x     = x->data;
  cl_uint cl_size = (cl_uint)y->size;

  cl_kernel kernel = op_opencl->kernels[0];

  /* Set kernel arguments */
  status = clSetKernelArg(kernel, 0, sizeof(cl_mem),  (void *)&cl_y);     AXB_ERRCHK(status);
  status = clSetKernelArg(kernel, 1, sizeof(cl_mem),  (void *)&cl_alpha); AXB_ERRCHK(status);
  status = clSetKernelArg(kernel, 2, sizeof(cl_mem),  (void *)&cl_x);     AXB_ERRCHK(status);
  status = clSetKernelArg(kernel, 3, sizeof(cl_uint), (void *)&cl_size);  AXB_ERRCHK(status);

  /* Run the kernel */
  size_t global_work_size = 256;
  size_t local_work_size = 256;
  status = clEnqueueNDRangeKernel(op_opencl->queue, kernel, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL); AXB_ERRCHK(status);
  status = clFinish(op_opencl->queue); AXB_ERRCHK(status);

  return 0;
}



static axbStatus_t axbOpBackendInit_OpenCL(axbHandle_t handle, axbOpBackend_t opencl_op_backend)
{
  axbMemBackend_t opencl_mem_backend;
  axbStatus_t status = axbMemBackendGetByName(handle, &opencl_mem_backend, "OpenCL");AXB_ERRCHK(status);

  opencl_op_backend->impl = malloc(sizeof(struct axbOpOpenCL_s));
  axbOpOpenCL_t op_opencl = (axbOpOpenCL_t)opencl_op_backend->impl;
  axbMemOpenCL_t mem_opencl = (axbMemOpenCL_t)opencl_mem_backend->impl;

  op_opencl->platform = mem_opencl->platform;
  op_opencl->context = mem_opencl->context;
  op_opencl->queue   = mem_opencl->queue;

  //
  // Create OpenCL programs:
  //
  op_opencl->programs_capacity = 1;
  op_opencl->programs = malloc(sizeof(cl_program) * op_opencl->programs_capacity);
  op_opencl->programs_size = 0;

  size_t kernel_len = strlen(my_compute_program);
  op_opencl->programs[op_opencl->programs_size] = clCreateProgramWithSource(op_opencl->context, 1, &my_compute_program, &kernel_len, &status); AXB_ERRCHK(status);
  op_opencl->programs_size += 1;

  status = clBuildProgram(op_opencl->programs[0], 1, mem_opencl->devices, NULL, NULL, NULL); AXB_ERRCHK(status);

  //
  // Create OpenCL kernels:
  // TODO: Think about reorganization. Kernels are attached to programs, not independent!
  //
  op_opencl->kernels_capacity = 1;
  op_opencl->kernels = malloc(sizeof(cl_kernel) * op_opencl->kernels_capacity);
  op_opencl->kernels_size = 0;

  op_opencl->kernels[op_opencl->kernels_size] = clCreateKernel(op_opencl->programs[0], "axpy", &status); AXB_ERRCHK(status);
  op_opencl->kernels_size += 1;

  return 0;
}

static axbStatus_t destroyOpenCLContext(void *impl)
{
  axbOpOpenCL_t opencl_context = (axbOpOpenCL_t)impl;
  for (size_t i=0; i<opencl_context->kernels_size; ++i) {
    cl_int status = clReleaseKernel(opencl_context->kernels[i]); AXB_ERRCHK(status);
  }
  free(opencl_context->kernels);
  for (size_t i=0; i<opencl_context->programs_size; ++i) {
    cl_int status = clReleaseProgram(opencl_context->programs[i]); AXB_ERRCHK(status);
  }
  free(opencl_context->programs);
  free(opencl_context);
  return 0;
}

axbStatus_t axbOpBackendRegister_OpenCL(axbHandle_t handle)
{
  axbOpBackend_t opencl_backend;
  axbStatus_t status = axbOpBackendCreate(&opencl_backend); AXB_ERRCHK(status);

  status = axbOpBackendInit_OpenCL(handle, opencl_backend); AXB_ERRCHK(status);

  // populate opencl_backend:
  status = axbOpBackendSetName(opencl_backend, "OpenCL"); AXB_ERRCHK(status);

  axbOperationID_t op_id = 0;
  status = axbOpBackendAddOperation(opencl_backend, "vec-axpy", (axbStatus_t (*)(void))op_axpy_opencl, opencl_backend->impl, &op_id); AXB_ERRCHK(status);

  status = axbOpBackendSetDestroy(opencl_backend, destroyOpenCLContext); AXB_ERRCHK(status);

  // push into enclosing context identified by handle:
  status = axbOpBackendRegister(handle, opencl_backend); AXB_ERRCHK(status);

  return 0;
}

#else

#include "libaxb/backend/opencl.h"

axbStatus_t axbOpBackendRegister_OpenCL(axbHandle_t handle)
{
  (void)handle;
  return 0;
}
#endif
