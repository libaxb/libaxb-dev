
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "libaxb.h"
#include "libaxb/backend.h"
#include "libaxb/general.h"
#include "libaxb/backend/op.h"

#include <cuda_runtime.h>

/////////////////////

__global__
static void kernel_vec_set_from_device(int n, double *x, const double *alpha)
{
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) x[i] = *alpha;
}

__global__
static void kernel_vec_set_from_host(int n, double *x, double alpha)
{
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) x[i] = alpha;
}


static axbStatus_t op_vec_set(axbVec_t x, axbScalar_t alpha, void *aux_data)
{
  (void)aux_data;

  double *d_x     = (double*)x->data;

  if (strcmp(alpha->memBackend->name, "host") != 0) {  // alpha on GPU
    double *d_alpha = (double*)alpha->data;
    kernel_vec_set_from_device<<<256, 256>>>((int)x->size, d_x, d_alpha);
  } else { // alpha on CPU
    double d_alpha = *((double*)alpha->data);
    kernel_vec_set_from_host<<<256, 256>>>((int)x->size, d_x, d_alpha);
  }

  return 0;
}

/////////////////////

__global__
static void kernel_vec_sqrtabs(int n, double *x)
{
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) x[i] = sqrt(fabs(x[i]));
}


static axbStatus_t op_vec_sqrtabs(axbVec_t x, void *aux_data)
{
  (void)aux_data;

  double *d_x     = (double*)x->data;
  kernel_vec_sqrtabs<<<256, 256>>>((int)x->size, d_x);

  return 0;
}

/////////////////////

static axbStatus_t op_vec_zero(axbVec_t x, void *aux_data)
{
  (void)aux_data;

  double *d_x     = (double*)x->data;
  kernel_vec_set_from_host<<<256, 256>>>((int)x->size, d_x, 0);

  return 0;
}

/////////////////////

__global__
static void kernel_vec_scale_from_device(int n, double *x, const double *alpha)
{
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) x[i] *= *alpha;
}

__global__
static void kernel_vec_scale_from_host(int n, double *x, double alpha)
{
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) x[i] *= alpha;
}


static axbStatus_t op_vec_scale(axbVec_t x, axbScalar_t alpha, void *aux_data)
{
  (void)aux_data;

  double *d_x     = (double*)x->data;

  if (strcmp(alpha->memBackend->name, "host") != 0) {  // alpha on GPU
    double *d_alpha = (double*)alpha->data;
    kernel_vec_scale_from_device<<<256, 256>>>((int)x->size, d_x, d_alpha);
  } else { // alpha on CPU
    double d_alpha = *((double*)alpha->data);
    kernel_vec_scale_from_host<<<256, 256>>>((int)x->size, d_x, d_alpha);
  }

  return 0;
}

//
// Reduction operations
//

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
__device__ double atomicAdd(double* address, double val) {
  unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif

/////////////////////

__global__
static void kernel_vec_sum(int n, double *x, double *alpha)
{
  __shared__ double reduction_buffer[256];
  double t = 0;
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) t += x[i];

  reduction_buffer[threadIdx.x] = t;

  // parallel reduction
  for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
  {
    __syncthreads();
    if (threadIdx.x < stride)
      reduction_buffer[threadIdx.x] += reduction_buffer[threadIdx.x+stride];
  }

  if (threadIdx.x == 0)
    atomicAdd(alpha, reduction_buffer[0]);
}


static axbStatus_t op_vec_sum(axbVec_t x, axbScalar_t alpha, void *aux_data)
{
  (void)aux_data;

  double *d_x     = (double*)x->data;
  double *d_alpha = (double*)alpha->data;
  kernel_vec_set_from_host<<<1, 1>>>((int)1, d_alpha, 0);
  kernel_vec_sum<<<256, 256>>>((int)x->size, d_x, d_alpha);

  return 0;
}

/////////////////////

__global__
static void kernel_vec_dot(int n, double *x, double *y, double *alpha)
{
  __shared__ double reduction_buffer[256];
  double t = 0;
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) t += x[i] * y[i];

  reduction_buffer[threadIdx.x] = t;

  // parallel reduction
  for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
  {
    __syncthreads();
    if (threadIdx.x < stride)
      reduction_buffer[threadIdx.x] += reduction_buffer[threadIdx.x+stride];
  }

  if (threadIdx.x == 0)
    atomicAdd(alpha, reduction_buffer[0]);
}


static axbStatus_t op_vec_dot(axbVec_t x, axbVec_t y, axbScalar_t alpha, void *aux_data)
{
  (void)aux_data;

  double *d_x     = (double*)x->data;
  double *d_y     = (double*)y->data;
  double *d_alpha = (double*)alpha->data;
  kernel_vec_set_from_host<<<1, 1>>>((int)1, d_alpha, 0);
  kernel_vec_dot<<<256, 256>>>((int)x->size, d_x, d_y, d_alpha);

  return 0;
}

/////////////////////

static axbStatus_t op_vec_tdot(axbVec_t x, axbVec_t y, axbScalar_t alpha, void *aux_data)
{
  return op_vec_dot(x, y, alpha, aux_data); // TODO: update for complex scalar types
}

/////////////////////

static axbStatus_t op_vec_mdot(axbVec_t x, size_t num_vecs, const axbVec_t *y, axbScalar_t *mdot, void *aux_data)
{
  (void)aux_data;

  // TODO: Replace by faster variant
  for (size_t i=0; i<num_vecs; ++i)
    op_vec_dot(x, y[i], mdot[i], aux_data);

  return 0;
}

/////////////////////

__global__
static void kernel_vec_norm1(int n, double *x, double *alpha)
{
  __shared__ double reduction_buffer[256];
  double t = 0;
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) t += fabs(x[i]);

  reduction_buffer[threadIdx.x] = t;

  // parallel reduction
  for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
  {
    __syncthreads();
    if (threadIdx.x < stride)
      reduction_buffer[threadIdx.x] += reduction_buffer[threadIdx.x+stride];
  }

  if (threadIdx.x == 0)
    atomicAdd(alpha, reduction_buffer[0]);
}


static axbStatus_t op_vec_norm1(axbVec_t x, axbScalar_t alpha, void *aux_data)
{
  (void)aux_data;

  double *d_x     = (double*)x->data;
  double *d_alpha = (double*)alpha->data;
  kernel_vec_set_from_host<<<1, 1>>>((int)1, d_alpha, 0);
  kernel_vec_norm1<<<256, 256>>>((int)x->size, d_x, d_alpha);

  return 0;
}

/////////////////////

__global__
static void kernel_vec_norm2(int n, double *x, double *alpha)
{
  __shared__ double reduction_buffer[256];
  double t = 0;
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) t += x[i] * x[i];

  reduction_buffer[threadIdx.x] = t;

  // parallel reduction
  for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
  {
    __syncthreads();
    if (threadIdx.x < stride)
      reduction_buffer[threadIdx.x] += reduction_buffer[threadIdx.x+stride];
  }

  if (threadIdx.x == 0)
    atomicAdd(alpha, reduction_buffer[0]);
}


static axbStatus_t op_vec_norm2(axbVec_t x, axbScalar_t alpha, void *aux_data)
{
  (void)aux_data;

  double *d_x     = (double*)x->data;
  double *d_alpha = (double*)alpha->data;
  kernel_vec_set_from_host<<<1, 1>>>((int)1, d_alpha, 0);
  kernel_vec_norm2<<<256, 256>>>((int)x->size, d_x, d_alpha);
  kernel_vec_sqrtabs<<<1, 1>>>((int)1, d_alpha);
  return 0;
}

/////////////////////

__global__
static void kernel_vec_norminf(int n, double *x, double *alpha)
{
  __shared__ double reduction_buffer[256];
  double t = 0;
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) t = max(t, fabs(x[i]));

  reduction_buffer[threadIdx.x] = t;

  // parallel reduction
  for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
  {
    __syncthreads();
    if (threadIdx.x < stride)
      reduction_buffer[threadIdx.x] = max(reduction_buffer[threadIdx.x], reduction_buffer[threadIdx.x+stride]);
  }

  if (threadIdx.x == 0)
    alpha[blockIdx.x] = reduction_buffer[0];
}


static axbStatus_t op_vec_norminf(axbVec_t x, axbScalar_t alpha, void *aux_data)
{
  (void)aux_data;

  double *tmp;
  cudaMalloc((void**)&tmp, sizeof(double) * 256);   // TODO: Avoid allocation in each call to op_vec_norminf

  double *d_x     = (double*)x->data;
  double *d_alpha = (double*)alpha->data;
  kernel_vec_set_from_host<<<1, 1>>>((int)1, d_alpha, 0);
  kernel_vec_norminf<<<256, 256>>>((int)x->size, d_x, tmp);
  kernel_vec_norminf<<<1, 256>>>((int)256, tmp, d_alpha);

  cudaFree(tmp);

  return 0;
}

/////////////////////

__global__
static void kernel_vec_dotnorm2(int n, double *s, double *t, double *dot_st, double *norm_t)
{
  __shared__ double reduction_buffer[256];
  double dot = 0;
  double norm = 0;
  *dot_st = 0;
  *norm_t = 0;

  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) {
    double val_t = t[i];
    dot += s[i] * val_t;
    norm += val_t * val_t;
  }

  //
  // first reduction for dot
  //
  reduction_buffer[threadIdx.x] = dot;

  // parallel reduction
  for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
  {
    __syncthreads();
    if (threadIdx.x < stride)
      reduction_buffer[threadIdx.x] += reduction_buffer[threadIdx.x+stride];
  }

  if (threadIdx.x == 0)
    atomicAdd(dot_st, reduction_buffer[0]);

  //
  // second reduction for norm
  //

  reduction_buffer[threadIdx.x] = norm;

  // parallel reduction
  for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
  {
    __syncthreads();
    if (threadIdx.x < stride)
      reduction_buffer[threadIdx.x] += reduction_buffer[threadIdx.x+stride];
  }

  if (threadIdx.x == 0)
    atomicAdd(norm_t, reduction_buffer[0]);
}


static axbStatus_t op_vec_dotnorm2(axbVec_t s, axbVec_t t, axbScalar_t dot, axbScalar_t norm, void *aux_data)
{
  (void)aux_data;

  double *d_s     = (double*)s->data;
  double *d_t     = (double*)t->data;
  double *d_dot   = (double*)dot->data;
  double *d_norm  = (double*)norm->data;
  kernel_vec_set_from_host<<<1, 1>>>((int)1, d_dot, 0);
  kernel_vec_set_from_host<<<1, 1>>>((int)1, d_norm, 0);
  kernel_vec_dotnorm2<<<256, 256>>>((int)s->size, d_s, d_t, d_dot, d_norm);
  kernel_vec_sqrtabs<<<1, 1>>>((int)1, d_norm);

  return 0;
}

/////////////////////

__global__
static void kernel_vec_max(int n, double *x, int *index, double *alpha)
{
  __shared__ double reduction_buffer_max[256];
  __shared__ int    reduction_buffer_idx[256];
  double t = x[0];
  int idx = 0;
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) {
    double xi = x[i];
    if (t < xi) {
      t = xi;
      idx = i;
    }
  }

  reduction_buffer_max[threadIdx.x] = t;
  reduction_buffer_idx[threadIdx.x] = idx;

  // parallel reduction
  for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
  {
    __syncthreads();
    if (threadIdx.x < stride) {
      if (reduction_buffer_max[threadIdx.x] < reduction_buffer_max[threadIdx.x+stride]) {
        reduction_buffer_max[threadIdx.x] = reduction_buffer_max[threadIdx.x+stride];
        reduction_buffer_idx[threadIdx.x] = reduction_buffer_idx[threadIdx.x+stride];
      }
    }
  }

  if (threadIdx.x == 0) {
    alpha[blockIdx.x] = reduction_buffer_max[0];
    index[blockIdx.x] = reduction_buffer_idx[0];
  }
}


static axbStatus_t op_vec_max(axbVec_t x, size_t *idx, axbScalar_t alpha, void *aux_data)
{
  (void)aux_data;

  //
  // TODO: Refactor this! The whole computation can be done without host<->device copies!
  //

  double *tmp;   cudaMalloc((void**)&tmp,   sizeof(double) * 256);   // TODO: Avoid allocation in each call to op_vec_max
  int    *index; cudaMalloc((void**)&index, sizeof(int)    * 256);   // TODO: Avoid allocation in each call to op_vec_max

  double *d_x     = (double*)x->data;
  double *d_alpha = (double*)alpha->data;
  kernel_vec_max<<<256, 256>>>((int)x->size, d_x, index, tmp);

  double host_val[256];
  cudaMemcpy(host_val, tmp,   256 * sizeof(double), cudaMemcpyDeviceToHost);

  int host_idx[256];
  cudaMemcpy(host_idx, index, 256 * sizeof(int), cudaMemcpyDeviceToHost);

  double val_max = host_val[0];
  int    idx_max = host_idx[0];
  for (size_t i=1; i<256; ++i) {
    if (val_max < host_val[i]) {
      val_max = host_val[i];
      idx_max = host_idx[i];
    }
  }
  *idx = idx_max;

  cudaMemcpy(d_alpha, (void*)&val_max, sizeof(double), cudaMemcpyHostToDevice);

  cudaFree(tmp);
  cudaFree(index);

  return 0;
}

/////////////////////

__global__
static void kernel_vec_min(int n, double *x, int *index, double *alpha)
{
  __shared__ double reduction_buffer_min[256];
  __shared__ int    reduction_buffer_idx[256];
  double t = x[0];
  int idx = 0;
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) {
    double xi = x[i];
    if (t > xi) {
      t = xi;
      idx = i;
    }
  }

  reduction_buffer_min[threadIdx.x] = t;
  reduction_buffer_idx[threadIdx.x] = idx;

  // parallel reduction
  for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
  {
    __syncthreads();
    if (threadIdx.x < stride) {
      if (reduction_buffer_min[threadIdx.x] > reduction_buffer_min[threadIdx.x+stride]) {
        reduction_buffer_min[threadIdx.x] = reduction_buffer_min[threadIdx.x+stride];
        reduction_buffer_idx[threadIdx.x] = reduction_buffer_idx[threadIdx.x+stride];
      }
    }
  }

  if (threadIdx.x == 0) {
    alpha[blockIdx.x] = reduction_buffer_min[0];
    index[blockIdx.x] = reduction_buffer_idx[0];
  }
}


static axbStatus_t op_vec_min(axbVec_t x, size_t *idx, axbScalar_t alpha, void *aux_data)
{
  (void)aux_data;

  //
  // TODO: Refactor this! The whole computation can be done without host<->device copies!
  //

  double *tmp;   cudaMalloc((void**)&tmp,   sizeof(double) * 256);   // TODO: Avoid allocation in each call to op_vec_max
  int    *index; cudaMalloc((void**)&index, sizeof(int)    * 256);   // TODO: Avoid allocation in each call to op_vec_max

  double *d_x     = (double*)x->data;
  double *d_alpha = (double*)alpha->data;
  kernel_vec_min<<<256, 256>>>((int)x->size, d_x, index, tmp);

  double host_val[256];
  cudaMemcpy(host_val, tmp,   256 * sizeof(double), cudaMemcpyDeviceToHost);

  int host_idx[256];
  cudaMemcpy(host_idx, index, 256 * sizeof(int), cudaMemcpyDeviceToHost);

  double val_min = host_val[0];
  int    idx_min = host_idx[0];
  for (size_t i=1; i<256; ++i) {
    if (val_min > host_val[i]) {
      val_min = host_val[i];
      idx_min = host_idx[i];
    }
  }
  *idx = idx_min;

  cudaMemcpy(d_alpha, (void*)&val_min, sizeof(double), cudaMemcpyHostToDevice);

  cudaFree(tmp);
  cudaFree(index);

  return 0;
}


//
// Vector-vector operations
//

__global__
static void kernel_vec_copy(int n, double *x, double *y)
{
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) y[i] = x[i];
}


static axbStatus_t op_vec_copy(axbVec_t x, axbVec_t y, void *aux_data)
{
  (void)aux_data;

  double *d_x     = (double*)x->data;
  double *d_y     = (double*)y->data;

  kernel_vec_copy<<<256, 256>>>((int)x->size, d_x, d_y);

  return 0;
}

/////////////////////

__global__
static void kernel_vec_swap(int n, double *x, double *y)
{
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) {
    double t = y[i];
    y[i] = x[i];
    x[i] = t;
  }
}


static axbStatus_t op_vec_swap(axbVec_t x, axbVec_t y, void *aux_data)
{
  (void)aux_data;

  double *d_x     = (double*)x->data;
  double *d_y     = (double*)y->data;

  kernel_vec_swap<<<256, 256>>>((int)x->size, d_x, d_y);

  return 0;
}

/////////////////////

__global__
static void kernel_vec_axpy(int n, double *y, const double *alpha, const double *x)
{
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) y[i] = *alpha * x[i] + y[i];
}


static axbStatus_t op_vec_axpy(axbVec_t y, axbScalar_t alpha, axbVec_t x, void *aux_data)
{
  (void)aux_data;

  double *d_y     = (double*)y->data;
  double *d_alpha = (double*)alpha->data;
  double *d_x     = (double*)x->data;

  kernel_vec_axpy<<<256, 256>>>((int)y->size, d_y, d_alpha, d_x);

  return 0;
}

/////////////////////

__global__
static void kernel_vec_aypx(int n, double *y, const double *alpha, const double *x)
{
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) y[i] = *alpha * y[i] + x[i];
}


static axbStatus_t op_vec_aypx(axbVec_t y, axbScalar_t alpha, axbVec_t x, void *aux_data)
{
  (void)aux_data;

  double *d_y     = (double*)y->data;
  double *d_alpha = (double*)alpha->data;
  double *d_x     = (double*)x->data;

  kernel_vec_aypx<<<256, 256>>>((int)y->size, d_y, d_alpha, d_x);

  return 0;
}

/////////////////////

__global__
static void kernel_vec_axpbypcz(int n, double *z, const double *alpha, const double *beta, const double *gamma, const double *x, const double *y)
{
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) z[i] = *alpha * x[i] + *beta * y[i] + *gamma * z[i];
}


static axbStatus_t op_vec_axpbypcz(axbVec_t z, axbScalar_t alpha, axbScalar_t beta, axbScalar_t gamma, axbVec_t x, axbVec_t y, void *aux_data)
{
  (void)aux_data;

  double *d_z     = (double*)z->data;
  double *d_alpha = (double*)alpha->data;
  double *d_beta  = (double*)beta->data;
  double *d_gamma = (double*)gamma->data;
  double *d_x     = (double*)x->data;
  double *d_y     = (double*)y->data;

  kernel_vec_axpbypcz<<<256, 256>>>((int)z->size, d_z, d_alpha, d_beta, d_gamma, d_x, d_y);

  return 0;
}

/////////////////////

__global__
static void kernel_vec_waxpy(int n, double *w, const double *alpha, const double *x, const double *y)
{
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) w[i] = *alpha * x[i] + y[i];
}


static axbStatus_t op_vec_waxpy(axbVec_t w, axbScalar_t alpha, axbVec_t x, axbVec_t y, void *aux_data)
{
  (void)aux_data;

  double *d_w     = (double*)w->data;
  double *d_alpha = (double*)alpha->data;
  double *d_x     = (double*)x->data;
  double *d_y     = (double*)y->data;

  kernel_vec_waxpy<<<256, 256>>>((int)w->size, d_w, d_alpha, d_x, d_y);

  return 0;
}

/////////////////////

static axbStatus_t op_vec_maxpy(axbVec_t y, size_t num_vecs, const axbScalar_t *alpha, const axbVec_t *x, void *aux_data) {

  // TODO: Be more efficient than this!
  for (size_t i=0; i<num_vecs; ++i)
    op_vec_axpy(y, alpha[i], x[i], aux_data);
  return 0;
}

/////////////////////

__global__
static void kernel_vec_pointwisemult(int n, double *w, const double *x, const double *y)
{
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) w[i] = x[i] * y[i];
}


static axbStatus_t op_vec_pointwisemult(axbVec_t w, axbVec_t x, axbVec_t y, void *aux_data)
{
  (void)aux_data;

  double *d_w     = (double*)w->data;
  double *d_x     = (double*)x->data;
  double *d_y     = (double*)y->data;

  kernel_vec_pointwisemult<<<256, 256>>>((int)w->size, d_w, d_x, d_y);

  return 0;
}

/////////////////////

__global__
static void kernel_vec_pointwisediv(int n, double *w, const double *x, const double *y)
{
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<n; i += gridDim.x * blockDim.x) w[i] = x[i] / y[i];
}


static axbStatus_t op_vec_pointwisediv(axbVec_t w, axbVec_t x, axbVec_t y, void *aux_data)
{
  (void)aux_data;

  double *d_w     = (double*)w->data;
  double *d_x     = (double*)x->data;
  double *d_y     = (double*)y->data;

  kernel_vec_pointwisediv<<<256, 256>>>((int)w->size, d_w, d_x, d_y);

  return 0;
}

/////////////////////



extern "C" axbStatus_t axbOpBackendRegister_CUDA(axbHandle_t handle)
{
  axbOpBackend_t cuda_backend;
  axbStatus_t status = axbOpBackendCreate(&cuda_backend); AXB_ERRCHK(status);

  // populate host_backend:
  status = axbOpBackendSetName(cuda_backend, "CUDA"); AXB_ERRCHK(status);

  axbOperationID_t op_id = 0;

#define AXB_ADD_OPERATION(OPNAME, ENUMCONSTANT)    status = axbOpBackendAddOperation(cuda_backend, #OPNAME,     (axbStatus_t (*)(void))OPNAME,     NULL, &op_id); AXB_ERRCHK(status); assert(op_id == ENUMCONSTANT && "Logic error: op_id != " #ENUMCONSTANT)

  // inplace operations
  AXB_ADD_OPERATION(op_vec_set,     AXB_OP_VEC_SET);
  AXB_ADD_OPERATION(op_vec_sqrtabs, AXB_OP_VEC_SQRTABS);
  AXB_ADD_OPERATION(op_vec_zero,    AXB_OP_VEC_ZERO);
  AXB_ADD_OPERATION(op_vec_scale,   AXB_OP_VEC_SCALE);

  // reduction operations
  AXB_ADD_OPERATION(op_vec_sum,      AXB_OP_VEC_SUM);
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
  AXB_ADD_OPERATION(op_vec_pointwisediv,  AXB_OP_VEC_POINTWISEDIV);

#undef AXB_ADD_OPERATION

  // push into enclosing context identified by handle:
  status = axbOpBackendRegister(handle, cuda_backend); AXB_ERRCHK(status);

  return 0;
}

