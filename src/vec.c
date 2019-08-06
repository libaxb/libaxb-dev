
#include <string.h>
#include <stdio.h>
#include "libaxb.h"
#include "libaxb/general.h"
#include "libaxb/backend/op.h"



#define AXB_VECTOR_INIT_CHECK(arg)     do { if (arg->init != 119346) { fprintf(stderr, "Vector not fully initialized!\n"); return 119346; } } while (0)

axbStatus_t axbVecCreateBegin(axbHandle_t handle, axbVec_t *vec)
{
  *vec = malloc(sizeof(struct axbVec_s));
  (*vec)->handle = handle;
  (*vec)->init = 119346;

  // set defaults:
  (*vec)->name = malloc(10);
  (*vec)->name[0] = 0;
  (*vec)->name_capacity = 10;
  (*vec)->data = NULL;

  (*vec)->datatype = AXB_REAL_DOUBLE;
  (*vec)->memBackend = handle->memBackends[0];
  (*vec)->opBackend  = handle->opBackends[0];

  return 0;
}

axbStatus_t axbVecSetSize(axbVec_t vec, size_t size)
{
  // TODO: Check that VecCreateEnd() has not been called yet!
  vec->size = size;
  return 0;
}
axbStatus_t axbVecGetSize(axbVec_t vec, size_t *size)
{
  *size = vec->size;
  return 0;
}

axbStatus_t axbVecSetDataType(axbVec_t vec, axbDataType_t datatype)
{
  // TODO: Check that VecCreateEnd() has not been called yet!
  vec->datatype = datatype;
  return 0;
}
axbStatus_t axbVecGetDataType(axbVec_t vec, axbDataType_t *datatype)
{
  *datatype = vec->datatype;
  return 0;
}

axbStatus_t axbVecSetMemBackend(axbVec_t vec, axbMemBackend_t backend)
{
  vec->memBackend = backend;
  return 0;
}
axbStatus_t axbVecGetMemBackend(axbVec_t vec, axbMemBackend_t *backend)
{
  *backend = vec->memBackend;
  return 0;
}

axbStatus_t axbVecSetOpBackend(axbVec_t vec, axbOpBackend_t backend)
{
  vec->opBackend = backend;
  return 0;
}
axbStatus_t axbVecGetOpBackend(axbVec_t vec, axbOpBackend_t *backend)
{
  *backend = vec->opBackend;
  return 0;
}


axbStatus_t axbVecCreateEnd(axbVec_t vec)
{
  if (vec->size > 0) {
    //printf("Initializing vector for backend %s\n", vec->memBackend->name);

    axbStatus_t status = vec->memBackend->op_malloc(&(vec->data), sizeof(double) * vec->size, vec->memBackend->impl);
    return status;
  }
  return 0;
}

axbStatus_t axbVecSetName(axbVec_t vec, const char *name)
{
  size_t len = strlen(name);
  if (len > vec->name_capacity) {
    free(vec->name);
    vec->name = malloc(len + 2);
  }
  for (size_t i=0; i<=len; ++i) vec->name[i] = name[i];
  return 0;
}
axbStatus_t axbVecGetName(axbVec_t vec, const char **name)
{
  *name = vec->name;
  return 0;
}

axbStatus_t axbVecSetValues(axbVec_t vec, void *values, axbDataType_t values_datatype)
{
  return axbMemBackendCopyIn(vec->memBackend, values, values_datatype, vec->data, vec->datatype, vec->size);
}
axbStatus_t axbVecGetValues(axbVec_t vec, void *values, axbDataType_t values_datatype)
{
  return axbMemBackendCopyOut(vec->memBackend, vec->data, vec->datatype, values, values_datatype, vec->size);
}

axbStatus_t axbVecDestroy(axbVec_t vec)
{
  AXB_VECTOR_INIT_CHECK(vec);
  vec->init += 1;

  vec->memBackend->op_free(vec->data, NULL);
  free(vec->name);
  free(vec);
  return 0;
}

//
// operations
//

// initialize vector:
axbStatus_t axbVecSet(axbVec_t x, axbScalar_t value)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_SET];
  axbStatus_t (*op)(axbVec_t, axbScalar_t, void*) = (axbStatus_t (*)(axbVec_t, axbScalar_t, void*)) op_desc.func;
  op(x, value, op_desc.func_data);
  return 0;
}

// in-place operations

axbStatus_t axbVecSqrtAbs(axbVec_t x)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_SQRTABS];
  axbStatus_t (*op)(axbVec_t, void*) = (axbStatus_t (*)(axbVec_t, void*)) op_desc.func;
  op(x, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecZero(axbVec_t x)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_ZERO];
  axbStatus_t (*op)(axbVec_t, void*) = (axbStatus_t (*)(axbVec_t, void*)) op_desc.func;
  op(x, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecScale(axbVec_t x, axbScalar_t alpha)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_SCALE];
  axbStatus_t (*op)(axbVec_t, axbScalar_t, void*) = (axbStatus_t (*)(axbVec_t, axbScalar_t, void*)) op_desc.func;
  op(x, alpha, op_desc.func_data);
  return 0;
}


// reduction operations:

axbStatus_t axbVecSum(axbVec_t x, axbScalar_t sum)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_SUM];
  axbStatus_t (*op)(axbVec_t, axbScalar_t, void*) = (axbStatus_t (*)(axbVec_t, axbScalar_t, void*)) op_desc.func;
  op(x, sum, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecDot(axbVec_t x, axbVec_t y, axbScalar_t dot)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_DOT];
  axbStatus_t (*op)(axbVec_t, axbVec_t, axbScalar_t, void*) = (axbStatus_t (*)(axbVec_t, axbVec_t, axbScalar_t, void*)) op_desc.func;
  op(x, y, dot, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecTDot(axbVec_t x, axbVec_t y, axbScalar_t tdot)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_TDOT];
  axbStatus_t (*op)(axbVec_t, axbVec_t, axbScalar_t, void*) = (axbStatus_t (*)(axbVec_t, axbVec_t, axbScalar_t, void*)) op_desc.func;
  op(x, y, tdot, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecMDot(axbVec_t x, size_t num_vecs, const axbVec_t *y, axbScalar_t *mdot)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_MDOT];
  axbStatus_t (*op)(axbVec_t, size_t, const axbVec_t *, axbScalar_t *, void*) = (axbStatus_t (*)(axbVec_t, size_t, const axbVec_t *, axbScalar_t *, void*)) op_desc.func;
  op(x, num_vecs, y, mdot, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecNorm1(axbVec_t x, axbScalar_t norm)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_NORM1];
  axbStatus_t (*op)(axbVec_t, axbScalar_t, void*) = (axbStatus_t (*)(axbVec_t, axbScalar_t, void*)) op_desc.func;
  op(x, norm, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecNorm2(axbVec_t x, axbScalar_t norm)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_NORM2];
  axbStatus_t (*op)(axbVec_t, axbScalar_t, void*) = (axbStatus_t (*)(axbVec_t, axbScalar_t, void*)) op_desc.func;
  op(x, norm, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecNormInf(axbVec_t x, axbScalar_t norm)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_NORMINF];
  axbStatus_t (*op)(axbVec_t, axbScalar_t, void*) = (axbStatus_t (*)(axbVec_t, axbScalar_t, void*)) op_desc.func;
  op(x, norm, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecDotNorm2(axbVec_t s, axbVec_t t, axbScalar_t dot_st, axbScalar_t norm_t)
{
  AXB_VECTOR_INIT_CHECK(s);
  AXB_VECTOR_INIT_CHECK(t);

  axbOpDescriptor_t op_desc = s->opBackend->op_table[AXB_OP_VEC_DOTNORM2];
  axbStatus_t (*op)(axbVec_t, axbVec_t, axbScalar_t, axbScalar_t, void*) = (axbStatus_t (*)(axbVec_t, axbVec_t, axbScalar_t, axbScalar_t, void*)) op_desc.func;
  op(s, t, dot_st, norm_t, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecMax(axbVec_t x, size_t *idx, axbScalar_t m)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_MAX];
  axbStatus_t (*op)(axbVec_t, size_t *, axbScalar_t, void*) = (axbStatus_t (*)(axbVec_t, size_t *, axbScalar_t, void*)) op_desc.func;
  op(x, idx, m, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecMin(axbVec_t x, size_t *idx, axbScalar_t m)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_MIN];
  axbStatus_t (*op)(axbVec_t, size_t *, axbScalar_t, void*) = (axbStatus_t (*)(axbVec_t, size_t *, axbScalar_t, void*)) op_desc.func;
  op(x, idx, m, op_desc.func_data);
  return 0;
}


// vector-vector operations:

axbStatus_t axbVecCopy(axbVec_t x, axbVec_t y)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);

  axbOpDescriptor_t op_desc = y->opBackend->op_table[AXB_OP_VEC_COPY];
  axbStatus_t (*op)(axbVec_t, axbVec_t, void*) = (axbStatus_t (*)(axbVec_t, axbVec_t, void*)) op_desc.func;
  op(x, y, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecSwap(axbVec_t x, axbVec_t y)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);

  axbOpDescriptor_t op_desc = y->opBackend->op_table[AXB_OP_VEC_SWAP];
  axbStatus_t (*op)(axbVec_t, axbVec_t, void*) = (axbStatus_t (*)(axbVec_t, axbVec_t, void*)) op_desc.func;
  op(x, y, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecAXPY(axbVec_t y, axbScalar_t alpha, axbVec_t x)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);

  axbOpDescriptor_t op_desc = y->opBackend->op_table[AXB_OP_VEC_AXPY];
  axbStatus_t (*op)(axbVec_t, axbScalar_t, axbVec_t, void*) = (axbStatus_t (*)(axbVec_t, axbScalar_t, axbVec_t, void*)) op_desc.func;
  op(y, alpha, x, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecAYPX(axbVec_t y, axbScalar_t alpha, axbVec_t x)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);

  axbOpDescriptor_t op_desc = y->opBackend->op_table[AXB_OP_VEC_AYPX];
  axbStatus_t (*op)(axbVec_t, axbScalar_t, axbVec_t, void*) = (axbStatus_t (*)(axbVec_t, axbScalar_t, axbVec_t, void*)) op_desc.func;
  op(y, alpha, x, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecAXPBYPCZ(axbVec_t z, axbScalar_t alpha, axbScalar_t beta, axbScalar_t gamma, axbVec_t x, axbVec_t y)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);
  AXB_VECTOR_INIT_CHECK(z);

  axbOpDescriptor_t op_desc = z->opBackend->op_table[AXB_OP_VEC_AXPBYPCZ];
  axbStatus_t (*op)(axbVec_t, axbScalar_t, axbScalar_t, axbScalar_t, axbVec_t, axbVec_t, void*) = (axbStatus_t (*)(axbVec_t, axbScalar_t, axbScalar_t, axbScalar_t, axbVec_t, axbVec_t, void*)) op_desc.func;
  op(z, alpha, beta, gamma, x, y, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecWAXPY(axbVec_t w, axbScalar_t alpha, axbVec_t x, axbVec_t y)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);
  AXB_VECTOR_INIT_CHECK(w);

  axbOpDescriptor_t op_desc = w->opBackend->op_table[AXB_OP_VEC_WAXPY];
  axbStatus_t (*op)(axbVec_t, axbScalar_t, axbVec_t, axbVec_t, void*) = (axbStatus_t (*)(axbVec_t, axbScalar_t, axbVec_t, axbVec_t, void*)) op_desc.func;
  op(w, alpha, x, y, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecMAXPY(axbVec_t y, size_t num_vecs, const axbScalar_t *alpha, const axbVec_t *x)
{
  AXB_VECTOR_INIT_CHECK(y);
  for (size_t i=0; i<num_vecs; ++i) AXB_VECTOR_INIT_CHECK(x[i]);

  axbOpDescriptor_t op_desc = y->opBackend->op_table[AXB_OP_VEC_MAXPY];
  axbStatus_t (*op)(axbVec_t, size_t, const axbScalar_t*, const axbVec_t*, void*) = (axbStatus_t (*)(axbVec_t, size_t, const axbScalar_t*, const axbVec_t*, void*)) op_desc.func;
  op(y, num_vecs, alpha, x, op_desc.func_data);
  return 0;
}



axbStatus_t axbVecPointwiseMult(axbVec_t w, axbVec_t x, axbVec_t y)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);
  AXB_VECTOR_INIT_CHECK(w);

  axbOpDescriptor_t op_desc = w->opBackend->op_table[AXB_OP_VEC_POINTWISEMULT];
  axbStatus_t (*op)(axbVec_t, axbVec_t, axbVec_t, void*) = (axbStatus_t (*)(axbVec_t, axbVec_t, axbVec_t, void*)) op_desc.func;
  op(w, x, y, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecPointwiseDivide(axbVec_t w, axbVec_t x, axbVec_t y)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);
  AXB_VECTOR_INIT_CHECK(w);

  axbOpDescriptor_t op_desc = w->opBackend->op_table[AXB_OP_VEC_POINTWISEDIV];
  axbStatus_t (*op)(axbVec_t, axbVec_t, axbVec_t, void*) = (axbStatus_t (*)(axbVec_t, axbVec_t, axbVec_t, void*)) op_desc.func;
  op(w, x, y, op_desc.func_data);
  return 0;
}


