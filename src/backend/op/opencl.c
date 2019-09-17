
#ifdef LIBAXB_ENABLE_OPENCL

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "libaxb.h"
#include "libaxb/backend.h"
#include "libaxb/general.h"
#include "libaxb/backend/op.h"
#include "libaxb/backend/opencl.h"

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
"__kernel void set_device(\n"
"          __global       double *y,\n"
"          __global const double *alpha, \n"
"          unsigned int size) \n"
"{ \n"
"  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0))\n"
"    y[i] = *alpha;\n"
"}\n"
"__kernel void set_host(\n"
"          __global double      *y,\n"
"                   double       alpha, \n"
"                   unsigned int size) \n"
"{ \n"
"  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0))\n"
"    y[i] = alpha;\n"
"}\n"
"__kernel void sqrtabs(\n"
"          __global double      *y,\n"
"                   unsigned int size) \n"
"{ \n"
"  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0))\n"
"    y[i] = sqrt(fabs(y[i]));\n"
"}\n"
"__kernel void scale_device(\n"
"          __global       double *y,\n"
"          __global const double *alpha, \n"
"          unsigned int size) \n"
"{ \n"
"  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0))\n"
"    y[i] *= *alpha;\n"
"}\n"
"__kernel void scale_host(\n"
"          __global double      *y,\n"
"                   double       alpha, \n"
"                   unsigned int size) \n"
"{ \n"
"  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0))\n"
"    y[i] *= alpha;\n"
"}\n"
"__kernel void axpy(\n"
"          __global       double *y,\n"
"          __global const double *alpha, \n"
"          __global const double *x,\n"
"          unsigned int size) \n"
"{ \n"
"  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0))\n"
"    y[i] += *alpha * x[i];\n"
"}\n"
"\n";


static axbStatus_t op_vec_set(struct axbVec_s *x, struct axbScalar_s *alpha, void *aux_data)
{
  cl_int status;
  struct axbOpOpenCL_s *op_opencl = (struct axbOpOpenCL_s*)aux_data;

  cl_mem cl_x     = x->data;
  cl_uint cl_size = (cl_uint)x->size;
  void *cl_alpha = alpha->data;
  size_t cl_alpha_size = sizeof(cl_mem);
  cl_kernel kernel = op_opencl->kernels[0];

  if (strcmp(alpha->memBackend->name, "host") != 0) {  // alpha on GPU
    kernel = op_opencl->kernels[1];
  } else {
    kernel = op_opencl->kernels[0];
    cl_alpha_size = sizeof(cl_double);
  }

  /* Set kernel arguments */
  status = clSetKernelArg(kernel, 0, sizeof(cl_mem),  (void *)&cl_x);     AXB_ERRCHK(status);
  status = clSetKernelArg(kernel, 1, cl_alpha_size,   (void *)&cl_alpha); AXB_ERRCHK(status);
  status = clSetKernelArg(kernel, 2, sizeof(cl_uint), (void *)&cl_size);  AXB_ERRCHK(status);

  /* Run the kernel */
  size_t global_work_size = 256;
  size_t local_work_size = 256;
  status = clEnqueueNDRangeKernel(op_opencl->queue, kernel, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL); AXB_ERRCHK(status);
  status = clFinish(op_opencl->queue); AXB_ERRCHK(status);

  return 0;
}

static axbStatus_t op_vec_sqrtabs(struct axbVec_s *x, void *aux_data)
{
  cl_int status;
  struct axbOpOpenCL_s *op_opencl = (struct axbOpOpenCL_s*)aux_data;

  cl_mem cl_x     = x->data;
  cl_uint cl_size = (cl_uint)x->size;
  cl_kernel kernel = op_opencl->kernels[2];

  /* Set kernel arguments */
  status = clSetKernelArg(kernel, 0, sizeof(cl_mem),  (void *)&cl_x);     AXB_ERRCHK(status);
  status = clSetKernelArg(kernel, 1, sizeof(cl_uint), (void *)&cl_size);  AXB_ERRCHK(status);

  /* Run the kernel */
  size_t global_work_size = 256;
  size_t local_work_size = 256;
  status = clEnqueueNDRangeKernel(op_opencl->queue, kernel, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL); AXB_ERRCHK(status);
  status = clFinish(op_opencl->queue); AXB_ERRCHK(status);

  return 0;
}

static axbStatus_t op_vec_zero(struct axbVec_s *x, void *aux_data)
{
  cl_int status;
  struct axbOpOpenCL_s *op_opencl = (struct axbOpOpenCL_s*)aux_data;

  cl_mem cl_x     = x->data;
  cl_uint cl_size = (cl_uint)x->size;
  double zero = 0;
  cl_kernel kernel = op_opencl->kernels[0];

  /* Set kernel arguments */
  status = clSetKernelArg(kernel, 0, sizeof(cl_mem),    (void *)&cl_x);    AXB_ERRCHK(status);
  status = clSetKernelArg(kernel, 1, sizeof(cl_double), (void *)&zero);    AXB_ERRCHK(status);
  status = clSetKernelArg(kernel, 2, sizeof(cl_uint),   (void *)&cl_size); AXB_ERRCHK(status);

  /* Run the kernel */
  size_t global_work_size = 256;
  size_t local_work_size = 256;
  status = clEnqueueNDRangeKernel(op_opencl->queue, kernel, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL); AXB_ERRCHK(status);
  status = clFinish(op_opencl->queue); AXB_ERRCHK(status);

  return 0;
}

static axbStatus_t op_vec_scale(struct axbVec_s *x, struct axbScalar_s *alpha, void *aux_data)
{
  cl_int status;
  struct axbOpOpenCL_s *op_opencl = (struct axbOpOpenCL_s*)aux_data;

  cl_mem cl_y     = x->data;
  cl_uint cl_size = (cl_uint)x->size;
  void *cl_alpha = alpha->data;
  size_t cl_alpha_size = sizeof(cl_mem);
  cl_kernel kernel;

  if (strcmp(alpha->memBackend->name, "host") != 0) {  // alpha on GPU
    kernel = op_opencl->kernels[4];
  } else {
    kernel = op_opencl->kernels[3];
    cl_alpha_size = sizeof(cl_double);
  }

  /* Set kernel arguments */
  status = clSetKernelArg(kernel, 0, sizeof(cl_mem),  (void *)&cl_y);     AXB_ERRCHK(status);
  status = clSetKernelArg(kernel, 1, cl_alpha_size,   (void *)&cl_alpha); AXB_ERRCHK(status);
  status = clSetKernelArg(kernel, 2, sizeof(cl_uint), (void *)&cl_size);  AXB_ERRCHK(status);

  /* Run the kernel */
  size_t global_work_size = 256;
  size_t local_work_size = 256;
  status = clEnqueueNDRangeKernel(op_opencl->queue, kernel, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL); AXB_ERRCHK(status);
  status = clFinish(op_opencl->queue); AXB_ERRCHK(status);

  return 0;
}

///////////////

static axbStatus_t op_vec_axpy(struct axbVec_s *y, struct axbScalar_s *alpha, struct axbVec_s *x, void *aux_data)
{
  cl_int status;
  struct axbOpOpenCL_s *op_opencl = (struct axbOpOpenCL_s*)aux_data;

  cl_mem cl_y     = y->data;
  cl_mem cl_alpha = alpha->data;
  cl_mem cl_x     = x->data;
  cl_uint cl_size = (cl_uint)y->size;

  cl_kernel kernel = op_opencl->kernels[2];

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



static axbStatus_t axbOpBackendInit_OpenCL(struct axbHandle_s *handle, struct axbOpBackend_s *opencl_op_backend)
{
  struct axbMemBackend_s *opencl_mem_backend;
  axbStatus_t status = axbMemBackendGetByName(handle, &opencl_mem_backend, "OpenCL");AXB_ERRCHK(status);

  opencl_op_backend->impl = malloc(sizeof(struct axbOpOpenCL_s));
  struct axbOpOpenCL_s *op_opencl = (struct axbOpOpenCL_s*)opencl_op_backend->impl;
  struct axbMemOpenCL_s *mem_opencl = (struct axbMemOpenCL_s*)opencl_mem_backend->impl;

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
  op_opencl->kernels_capacity = 32;
  op_opencl->kernels = malloc(sizeof(cl_kernel) * op_opencl->kernels_capacity);
  op_opencl->kernels_size = 0;

  op_opencl->kernels[op_opencl->kernels_size] = clCreateKernel(op_opencl->programs[0], "set_host", &status); AXB_ERRCHK(status);      // 0
  op_opencl->kernels_size += 1;
  op_opencl->kernels[op_opencl->kernels_size] = clCreateKernel(op_opencl->programs[0], "set_device", &status); AXB_ERRCHK(status);
  op_opencl->kernels_size += 1;
  op_opencl->kernels[op_opencl->kernels_size] = clCreateKernel(op_opencl->programs[0], "sqrtabs", &status); AXB_ERRCHK(status);
  op_opencl->kernels_size += 1;
  op_opencl->kernels[op_opencl->kernels_size] = clCreateKernel(op_opencl->programs[0], "scale_host", &status); AXB_ERRCHK(status);
  op_opencl->kernels_size += 1;
  op_opencl->kernels[op_opencl->kernels_size] = clCreateKernel(op_opencl->programs[0], "scale_device", &status); AXB_ERRCHK(status);
  op_opencl->kernels_size += 1;

  op_opencl->kernels[op_opencl->kernels_size] = clCreateKernel(op_opencl->programs[0], "axpy", &status); AXB_ERRCHK(status);
  op_opencl->kernels_size += 1;

  return 0;
}

static axbStatus_t destroyOpenCLContext(void *impl)
{
  struct axbOpOpenCL_s *opencl_context = (struct axbOpOpenCL_s*)impl;
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

axbStatus_t axbOpBackendRegister_OpenCL(struct axbHandle_s *handle)
{
  struct axbOpBackend_s *opencl_backend;
  axbStatus_t status = axbOpBackendCreate(&opencl_backend); AXB_ERRCHK(status);

  status = axbOpBackendInit_OpenCL(handle, opencl_backend); AXB_ERRCHK(status);

  // populate opencl_backend:
  status = axbOpBackendSetName(opencl_backend, "OpenCL"); AXB_ERRCHK(status);

  axbOperationID_t op_id = 0;

#define AXB_ADD_OPERATION(OPNAME, ENUMCONSTANT)    status = axbOpBackendAddOperation(opencl_backend, #OPNAME, (axbStatus_t (*)(void))OPNAME, opencl_backend->impl, &op_id); AXB_ERRCHK(status); assert(op_id == ENUMCONSTANT && "Logic error: op_id != " #ENUMCONSTANT)

  // inplace operations
  AXB_ADD_OPERATION(op_vec_set,     AXB_OP_VEC_SET);
  AXB_ADD_OPERATION(op_vec_sqrtabs, AXB_OP_VEC_SQRTABS);
  AXB_ADD_OPERATION(op_vec_zero,    AXB_OP_VEC_ZERO);
  AXB_ADD_OPERATION(op_vec_scale,   AXB_OP_VEC_SCALE);

  // reduction operations
  /*AXB_ADD_OPERATION(op_vec_sum,      AXB_OP_VEC_SUM);
  AXB_ADD_OPERATION(op_vec_dot,      AXB_OP_VEC_DOT);
  AXB_ADD_OPERATION(op_vec_tdot,     AXB_OP_VEC_TDOT);
  AXB_ADD_OPERATION(op_vec_mdot,     AXB_OP_VEC_MDOT);
  AXB_ADD_OPERATION(op_vec_norm1,    AXB_OP_VEC_NORM1);
  AXB_ADD_OPERATION(op_vec_norm2,    AXB_OP_VEC_NORM2);
  AXB_ADD_OPERATION(op_vec_norminf,  AXB_OP_VEC_NORMINF);
  AXB_ADD_OPERATION(op_vec_dotnorm2, AXB_OP_VEC_DOTNORM2);
  AXB_ADD_OPERATION(op_vec_max,      AXB_OP_VEC_MAX);
  AXB_ADD_OPERATION(op_vec_min,      AXB_OP_VEC_MIN);

  // vector-vector operations
  AXB_ADD_OPERATION(op_vec_copy,          AXB_OP_VEC_COPY);
  AXB_ADD_OPERATION(op_vec_swap,          AXB_OP_VEC_SWAP);
  AXB_ADD_OPERATION(op_vec_axpy,          AXB_OP_VEC_AXPY);
  AXB_ADD_OPERATION(op_vec_aypx,          AXB_OP_VEC_AYPX);
  AXB_ADD_OPERATION(op_vec_axpbypcz,      AXB_OP_VEC_AXPBYPCZ);
  AXB_ADD_OPERATION(op_vec_waxpy,         AXB_OP_VEC_WAXPY);
  AXB_ADD_OPERATION(op_vec_maxpy,         AXB_OP_VEC_MAXPY);
  AXB_ADD_OPERATION(op_vec_pointwisemult, AXB_OP_VEC_POINTWISEMULT);
  AXB_ADD_OPERATION(op_vec_pointwisediv,  AXB_OP_VEC_POINTWISEDIV);*/

#undef AXB_ADD_OPERATION

  status = axbOpBackendSetDestroy(opencl_backend, destroyOpenCLContext); AXB_ERRCHK(status);

  // push into enclosing context identified by handle:
  status = axbOpBackendRegister(handle, opencl_backend); AXB_ERRCHK(status);

  return 0;
}

#else

#include "libaxb/backend/opencl.h"

axbStatus_t axbOpBackendRegister_OpenCL(struct axbHandle_s *handle)
{
  (void)handle;
  return 0;
}
#endif
