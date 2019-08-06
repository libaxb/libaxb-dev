
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "libaxb.h"
#include "libaxb/backend.h"
#include "libaxb/general.h"
#include "libaxb/backend/op.h"

/*
 * General note:
 *  All op_* routines are worker routines and need to be reimplemented for each op-backend.
 *  Backend-agnostic checks like correct data types or compatible vector sizes are assumed to be checked at a higher level (->backend-agnostic), rather than being duplicated in each worker routine.
 */

static axbStatus_t op_vec_set(axbVec_t x, axbScalar_t alpha, void *aux_data)
{
  (void)aux_data;

  double *d_x     = x->data;
  double *d_alpha = (double*) alpha->data;

  for (size_t i=0; i<x->size; ++i) d_x[i] = *d_alpha;

  return 0;
}

static axbStatus_t op_vec_sqrtabs(axbVec_t x, void *aux_data)
{
  (void)aux_data;

  double *d_x     = x->data;

  for (size_t i=0; i<x->size; ++i) d_x[i] = sqrt(fabs(d_x[i]));

  return 0;
}

static axbStatus_t op_vec_zero(axbVec_t x, void *aux_data)
{
  (void)aux_data;

  double *d_x     = x->data;

  for (size_t i=0; i<x->size; ++i) d_x[i] = 0;

  return 0;
}

static axbStatus_t op_vec_scale(axbVec_t x, axbScalar_t alpha, void *aux_data)
{
  (void)aux_data;

  double *d_x     = x->data;
  double *d_alpha = (double*) alpha->data;

  for (size_t i=0; i<x->size; ++i) d_x[i] *= *d_alpha;

  return 0;
}

//////// reductions

static axbStatus_t op_vec_sum(axbVec_t x, axbScalar_t sum, void *aux_data)
{
  assert(sum->datatype == AXB_REAL_DOUBLE);
  (void)aux_data;

  double *d_x   = x->data;
  double *d_sum = (double*) sum->data;

  double s = 0;
  for (size_t i=0; i<x->size; ++i) s += d_x[i];
  *d_sum = s;

  return 0;
}

static axbStatus_t op_vec_dot(axbVec_t x, axbVec_t y, axbScalar_t dot, void *aux_data)
{
  (void)aux_data;

  double *d_x     = x->data;
  double *d_y     = y->data;
  double *d_dot = (double*) dot->data;

  double d = 0;
  for (size_t i=0; i<y->size; ++i) d += d_x[i] * d_y[i];
  *d_dot = d;

  return 0;
}

static axbStatus_t op_vec_tdot(axbVec_t x, axbVec_t y, axbScalar_t tdot, void *aux_data)
{
  (void)aux_data;

  double *d_x     = x->data;
  double *d_y     = y->data;
  double *d_tdot = (double*) tdot->data;

  double t = 0;
  for (size_t i=0; i<x->size; ++i) t += d_x[i] * d_y[i];
  *d_tdot = t;

  return 0;
}

static axbStatus_t op_vec_mdot(axbVec_t x, size_t num_vecs, const axbVec_t *y, axbScalar_t *mdot, void *aux_data)
{
  (void)aux_data;

  double *d_x    = x->data;
  double val_x = 0;

  // auxiliary variables
  double m[4] = {0};

  while (num_vecs) {
    m[0] = m[1] = m[2] = m[3] = 0;

    // 'unrolling' for up to four right hand sides y simultaneously.
    // This leverages most of the performance benefits by reducing memory transfers.
    // Further optimizations by unrolling up to 8 vectors are possible.
    // Unrolling to 9+ vectors results in diminishing returns (<10 percent gain)
    switch (num_vecs) {
    case 3:
      for (size_t i=0; i<x->size; ++i) {
        val_x = d_x[i];
        m[0] += val_x * ((double*)y[0])[i];
        m[1] += val_x * ((double*)y[1])[i];
        m[2] += val_x * ((double*)y[2])[i];
      }
      *((double*)mdot[0]->data) = m[0];
      *((double*)mdot[1]->data) = m[1];
      *((double*)mdot[2]->data) = m[2];

      num_vecs = 0;
      break;
    case 2:
      for (size_t i=0; i<x->size; ++i) {
        val_x = d_x[i];
        m[0] += val_x * ((double*)y[0])[i];
        m[1] += val_x * ((double*)y[1])[i];
      }
      *((double*)mdot[0]->data) = m[0];
      *((double*)mdot[1]->data) = m[1];

      num_vecs = 0;
      break;
    case 1:
      for (size_t i=0; i<x->size; ++i) {
        m[0] += d_x[i] * ((double*)y[0])[i];
      }
      *((double*)mdot[0]->data) = m[0];

      num_vecs = 0;
      break;
    default:
      for (size_t i=0; i<x->size; ++i) {
        val_x = d_x[i];
        m[0] += val_x * ((double*)y[0])[i];
        m[1] += val_x * ((double*)y[1])[i];
        m[2] += val_x * ((double*)y[2])[i];
        m[3] += val_x * ((double*)y[3])[i];
      }
      *((double*)mdot[0]->data) = m[0];
      *((double*)mdot[1]->data) = m[1];
      *((double*)mdot[2]->data) = m[2];
      *((double*)mdot[3]->data) = m[3];

      y += 4;
      mdot += 4;
      num_vecs -= 4;
      break;
    }
  } // while

  return 0;
}

static axbStatus_t op_vec_norm1(axbVec_t x, axbScalar_t norm, void *aux_data)
{
  (void)aux_data;

  double *d_x     = x->data;
  double *d_norm = (double*) norm->data;

  double n = 0;
  for (size_t i=0; i<x->size; ++i) n += fabs(d_x[i]);
  *d_norm = n;

  return 0;
}

static axbStatus_t op_vec_norm2(axbVec_t x, axbScalar_t norm, void *aux_data)
{
  (void)aux_data;

  double *d_x     = x->data;
  double *d_norm = (double*) norm->data;

  double n = 0;
  for (size_t i=0; i<x->size; ++i) n += d_x[i] * d_x[i];
  *d_norm = sqrt(n);

  return 0;
}

static axbStatus_t op_vec_norminf(axbVec_t x, axbScalar_t norm, void *aux_data)
{
  (void)aux_data;

  double *d_x     = x->data;
  double *d_norm = (double*) norm->data;

  double n = 0;
  for (size_t i=0; i<x->size; ++i) {
    if (fabs(d_x[i]) > n) n = fabs(d_x[i]);
  }
  *d_norm = n;

  return 0;
}

static axbStatus_t op_vec_dotnorm2(axbVec_t s, axbVec_t t, axbScalar_t dot_st, axbScalar_t norm_t, void *aux_data)
{
  (void)aux_data;

  double *d_s    = s->data;
  double *d_t    = t->data;
  double *d_dot  = (double*) dot_st->data;
  double *d_norm = (double*) norm_t->data;

  double d = 0;
  double n = 0;
  for (size_t i=0; i<s->size; ++i) {
    d += d_s[i] * d_t[i];
    n += d_t[i] * d_t[i];
  }
  *d_dot  = d;
  *d_norm = sqrt(n);

  return 0;
}

static axbStatus_t op_vec_max(axbVec_t x, size_t *idx, axbScalar_t m, void *aux_data)
{
  (void)aux_data;

  double *d_x   = x->data;
  double *d_max = (double*) m->data;

  double cur_max = d_x[0];
  *idx = 0;
  for (size_t i=1; i<x->size; ++i) {
    if (d_x[i] > cur_max) {
      cur_max = d_x[i];
      *idx = i;
    }
  }
  *d_max = cur_max;

  return 0;
}

static axbStatus_t op_vec_min(axbVec_t x, size_t *idx, axbScalar_t m, void *aux_data)
{
  (void)aux_data;

  double *d_x   = x->data;
  double *d_min = (double*) m->data;

  double cur_min = d_x[0];
  *idx = 0;
  for (size_t i=1; i<x->size; ++i) {
    if (d_x[i] < cur_min) {
      cur_min = d_x[i];
      *idx = i;
    }
  }
  *d_min = cur_min;

  return 0;
}

/////////////// vector-vector

static axbStatus_t op_vec_copy(axbVec_t x, axbVec_t y, void *aux_data)
{
  (void)aux_data;

  double *d_x     = x->data;
  double *d_y     = y->data;

  for (size_t i=0; i<x->size; ++i) d_y[i] = d_x[i];

  return 0;
}

static axbStatus_t op_vec_swap(axbVec_t x, axbVec_t y, void *aux_data)
{
  (void)aux_data;

  double *d_x     = x->data;
  double *d_y     = y->data;

  for (size_t i=0; i<y->size; ++i) {
    double t = d_x[i];
    d_x[i] = d_y[i];
    d_y[i] = t;
  }

  return 0;
}

static axbStatus_t op_vec_axpy(axbVec_t y, axbScalar_t alpha, axbVec_t x, void *aux_data)
{
  (void)aux_data;

  double *d_y     = y->data;
  double *d_alpha = (double*) alpha->data;
  double *d_x     = x->data;

  for (size_t i=0; i<y->size; ++i) d_y[i] += *d_alpha * d_x[i];

  return 0;
}

static axbStatus_t op_vec_aypx(axbVec_t y, axbScalar_t alpha, axbVec_t x, void *aux_data)
{
  (void)aux_data;

  double *d_y     = y->data;
  double *d_alpha = (double*) alpha->data;
  double *d_x     = x->data;

  for (size_t i=0; i<y->size; ++i) d_y[i] = *d_alpha * d_y[i] + d_x[i];

  return 0;
}

static axbStatus_t op_vec_axpbypcz(axbVec_t z, axbScalar_t alpha, axbScalar_t beta, axbScalar_t gamma, axbVec_t x, axbVec_t y, void *aux_data)
{
  (void)aux_data;

  double *d_z     = z->data;
  double *d_alpha = (double*) alpha->data;
  double *d_beta  = (double*) beta->data;
  double *d_gamma = (double*) gamma->data;
  double *d_x     = x->data;
  double *d_y     = y->data;

  for (size_t i=0; i<z->size; ++i) d_z[i] = *d_alpha * d_x[i] + *d_beta * d_y[i] + *d_gamma * d_z[i];

  return 0;
}

static axbStatus_t op_vec_waxpy(axbVec_t w, axbScalar_t alpha, axbVec_t x, axbVec_t y, void *aux_data)
{
  (void)aux_data;

  double *d_w     = w->data;
  double *d_alpha = (double*) alpha->data;
  double *d_x     = x->data;
  double *d_y     = y->data;

  for (size_t i=0; i<y->size; ++i) d_w[i] = *d_alpha * d_x[i] + d_y[i];

  return 0;
}

static axbStatus_t op_vec_maxpy(axbVec_t y, size_t num_vecs, const axbScalar_t *alpha, const axbVec_t *x, void *aux_data)
{
  (void)aux_data;

  double *d_y     = y->data;

  for (size_t i=0; i<y->size; ++i) {
    double ax = 0;
    for (size_t j=0; j<num_vecs; ++j) {
      double *d_alpha = (double*) alpha[j]->data;
      double *d_x     = x[j]->data;

      ax += *d_alpha * d_x[i];
    }
    d_y[i] += ax;
  }

  return 0;
}

static axbStatus_t op_vec_pointwisemult(axbVec_t w, axbVec_t x, axbVec_t y, void *aux_data)
{
  (void)aux_data;

  double *d_w     = w->data;
  double *d_x     = x->data;
  double *d_y     = y->data;

  for (size_t i=0; i<w->size; ++i) d_w[i] = d_x[i] * d_y[i];

  return 0;
}

static axbStatus_t op_vec_pointwisediv(axbVec_t w, axbVec_t x, axbVec_t y, void *aux_data)
{
  (void)aux_data;

  double *d_w     = w->data;
  double *d_x     = x->data;
  double *d_y     = y->data;

  for (size_t i=0; i<w->size; ++i) d_w[i] = d_x[i] / d_y[i];

  return 0;
}


//
// Registration of routines
//

axbStatus_t axbOpBackendRegister_Host(axbHandle_t handle)
{
  axbOpBackend_t host_backend;
  axbStatus_t status = axbOpBackendCreate(&host_backend); AXB_ERRCHK(status);

  // populate host_backend:
  status = axbOpBackendSetName(host_backend, "host"); AXB_ERRCHK(status);

  axbOperationID_t op_id = 0;

#define AXB_ADD_OPERATION(OPNAME, ENUMCONSTANT)    status = axbOpBackendAddOperation(host_backend, #OPNAME,     (axbStatus_t (*)(void))OPNAME,     NULL, &op_id); AXB_ERRCHK(status); assert(op_id == ENUMCONSTANT && "Logic error: op_id != " #ENUMCONSTANT)

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
  status = axbOpBackendRegister(handle, host_backend); AXB_ERRCHK(status);

  return 0;
}

