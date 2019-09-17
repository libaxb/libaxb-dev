
#include <string.h>
#include <stdio.h>
#include "libaxb.h"
#include "libaxb/general.h"
#include "libaxb/backend/op.h"



#define AXB_VECTOR_INIT_CHECK(arg)     do { if (arg->init != 119346) { fprintf(stderr, "Vector not fully initialized!\n"); return 119346; } } while (0)

axbStatus_t axbVecCreateBegin(struct axbHandle_s *handle, struct axbVec_s **vec)
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

axbStatus_t axbVecSetSize(struct axbVec_s *vec, size_t size)
{
  // TODO: Check that VecCreateEnd() has not been called yet!
  vec->size = size;
  return 0;
}
axbStatus_t axbVecGetSize(const struct axbVec_s *vec, size_t *size)
{
  *size = vec->size;
  return 0;
}

axbStatus_t axbVecSetDataType(struct axbVec_s *vec, axbDataType_t datatype)
{
  // TODO: Check that VecCreateEnd() has not been called yet!
  vec->datatype = datatype;
  return 0;
}
axbStatus_t axbVecGetDataType(const struct axbVec_s *vec, axbDataType_t *datatype)
{
  *datatype = vec->datatype;
  return 0;
}

axbStatus_t axbVecSetMemBackend(struct axbVec_s *vec, struct axbMemBackend_s *backend)
{
  vec->memBackend = backend;
  return 0;
}
axbStatus_t axbVecGetMemBackend(const struct axbVec_s *vec, struct axbMemBackend_s **backend)
{
  *backend = vec->memBackend;
  return 0;
}

axbStatus_t axbVecSetOpBackend(struct axbVec_s *vec, struct axbOpBackend_s *backend)
{
  vec->opBackend = backend;
  return 0;
}
axbStatus_t axbVecGetOpBackend(const struct axbVec_s *vec, struct axbOpBackend_s **backend)
{
  *backend = vec->opBackend;
  return 0;
}


axbStatus_t axbVecCreateEnd(struct axbVec_s *vec)
{
  if (vec->size > 0) {
    //printf("Initializing vector for backend %s\n", vec->memBackend->name);

    axbStatus_t status = vec->memBackend->op_malloc(&(vec->data), sizeof(double) * vec->size, vec->memBackend->impl);
    return status;
  }
  return 0;
}

axbStatus_t axbVecSetName(struct axbVec_s *vec, const char *name)
{
  size_t len = strlen(name);
  if (len > vec->name_capacity) {
    free(vec->name);
    vec->name = malloc(len + 2);
  }
  for (size_t i=0; i<=len; ++i) vec->name[i] = name[i];
  return 0;
}
axbStatus_t axbVecGetName(const struct axbVec_s *vec, const char **name)
{
  *name = vec->name;
  return 0;
}

axbStatus_t axbVecSetValues(struct axbVec_s *vec, void *values, axbDataType_t values_datatype)
{
  return axbMemBackendCopyIn(vec->memBackend, values, values_datatype, vec->data, vec->datatype, vec->size);
}
axbStatus_t axbVecGetValues(const struct axbVec_s *vec, void *values, axbDataType_t values_datatype)
{
  return axbMemBackendCopyOut(vec->memBackend, vec->data, vec->datatype, values, values_datatype, vec->size);
}

axbStatus_t axbVecDestroy(struct axbVec_s *vec)
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
axbStatus_t axbVecSet(struct axbVec_s *x, const struct axbScalar_s *value)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_SET];
  axbStatus_t (*op)(struct axbVec_s *, const struct axbScalar_s *, void*) = (axbStatus_t (*)(struct axbVec_s *, const struct axbScalar_s *, void*)) op_desc.func;
  op(x, value, op_desc.func_data);
  return 0;
}

// in-place operations

axbStatus_t axbVecSqrtAbs(struct axbVec_s *x)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_SQRTABS];
  axbStatus_t (*op)(struct axbVec_s *, void*) = (axbStatus_t (*)(struct axbVec_s *, void*)) op_desc.func;
  op(x, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecZero(struct axbVec_s *x)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_ZERO];
  axbStatus_t (*op)(struct axbVec_s *, void*) = (axbStatus_t (*)(struct axbVec_s *, void*)) op_desc.func;
  op(x, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecScale(struct axbVec_s *x, const struct axbScalar_s *alpha)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_SCALE];
  axbStatus_t (*op)(struct axbVec_s *, const struct axbScalar_s *, void*) = (axbStatus_t (*)(struct axbVec_s *, const struct axbScalar_s *, void*)) op_desc.func;
  op(x, alpha, op_desc.func_data);
  return 0;
}


// reduction operations:

axbStatus_t axbVecSum(const struct axbVec_s *x, struct axbScalar_s *sum)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_SUM];
  axbStatus_t (*op)(const struct axbVec_s *, struct axbScalar_s *, void*) = (axbStatus_t (*)(const struct axbVec_s *, struct axbScalar_s *, void*)) op_desc.func;
  op(x, sum, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecDot(const struct axbVec_s *x, const struct axbVec_s *y, struct axbScalar_s *dot)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_DOT];
  axbStatus_t (*op)(const struct axbVec_s *, const struct axbVec_s *, struct axbScalar_s *, void*) = (axbStatus_t (*)(const struct axbVec_s *, const struct axbVec_s *, struct axbScalar_s *, void*)) op_desc.func;
  op(x, y, dot, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecTDot(const struct axbVec_s *x, const struct axbVec_s *y, struct axbScalar_s *tdot)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_TDOT];
  axbStatus_t (*op)(const struct axbVec_s *, const struct axbVec_s *, struct axbScalar_s *, void*) = (axbStatus_t (*)(const struct axbVec_s *, const struct axbVec_s *, struct axbScalar_s *, void*)) op_desc.func;
  op(x, y, tdot, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecMDot(const struct axbVec_s *x, size_t num_vecs, const struct axbVec_s **y, struct axbScalar_s **mdot)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_MDOT];
  axbStatus_t (*op)(const struct axbVec_s *, size_t, const struct axbVec_s **, struct axbScalar_s **, void*) = (axbStatus_t (*)(const struct axbVec_s *, size_t, const struct axbVec_s **, struct axbScalar_s **, void*)) op_desc.func;
  op(x, num_vecs, y, mdot, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecNorm1(const struct axbVec_s *x, struct axbScalar_s *norm)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_NORM1];
  axbStatus_t (*op)(const struct axbVec_s *, struct axbScalar_s *, void*) = (axbStatus_t (*)(const struct axbVec_s *, struct axbScalar_s *, void*)) op_desc.func;
  op(x, norm, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecNorm2(const struct axbVec_s *x, struct axbScalar_s *norm)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_NORM2];
  axbStatus_t (*op)(const struct axbVec_s *, struct axbScalar_s *, void*) = (axbStatus_t (*)(const struct axbVec_s *, struct axbScalar_s *, void*)) op_desc.func;
  op(x, norm, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecNormInf(const struct axbVec_s *x, struct axbScalar_s *norm)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_NORMINF];
  axbStatus_t (*op)(const struct axbVec_s *, struct axbScalar_s *, void*) = (axbStatus_t (*)(const struct axbVec_s *, struct axbScalar_s *, void*)) op_desc.func;
  op(x, norm, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecDotNorm2(const struct axbVec_s *s, const struct axbVec_s *t, struct axbScalar_s *dot_st, struct axbScalar_s *norm_t)
{
  AXB_VECTOR_INIT_CHECK(s);
  AXB_VECTOR_INIT_CHECK(t);

  axbOpDescriptor_t op_desc = s->opBackend->op_table[AXB_OP_VEC_DOTNORM2];
  axbStatus_t (*op)(const struct axbVec_s *, const struct axbVec_s *, struct axbScalar_s *, struct axbScalar_s *, void*) = (axbStatus_t (*)(const struct axbVec_s *, const struct axbVec_s *, struct axbScalar_s *, struct axbScalar_s *, void*)) op_desc.func;
  op(s, t, dot_st, norm_t, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecMax(const struct axbVec_s *x, size_t *idx, struct axbScalar_s *m)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_MAX];
  axbStatus_t (*op)(const struct axbVec_s *, size_t *, struct axbScalar_s *, void*) = (axbStatus_t (*)(const struct axbVec_s *, size_t *, struct axbScalar_s *, void*)) op_desc.func;
  op(x, idx, m, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecMin(const struct axbVec_s *x, size_t *idx, struct axbScalar_s *m)
{
  AXB_VECTOR_INIT_CHECK(x);

  axbOpDescriptor_t op_desc = x->opBackend->op_table[AXB_OP_VEC_MIN];
  axbStatus_t (*op)(const struct axbVec_s *, size_t *, struct axbScalar_s *, void*) = (axbStatus_t (*)(const struct axbVec_s *, size_t *, struct axbScalar_s *, void*)) op_desc.func;
  op(x, idx, m, op_desc.func_data);
  return 0;
}


// vector-vector operations:

axbStatus_t axbVecCopy(const struct axbVec_s *x, struct axbVec_s *y)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);

  axbOpDescriptor_t op_desc = y->opBackend->op_table[AXB_OP_VEC_COPY];
  axbStatus_t (*op)(const struct axbVec_s *, struct axbVec_s *, void*) = (axbStatus_t (*)(const struct axbVec_s *, struct axbVec_s *, void*)) op_desc.func;
  op(x, y, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecSwap(struct axbVec_s *x, struct axbVec_s *y)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);

  axbOpDescriptor_t op_desc = y->opBackend->op_table[AXB_OP_VEC_SWAP];
  axbStatus_t (*op)(struct axbVec_s *, struct axbVec_s *, void*) = (axbStatus_t (*)(struct axbVec_s *, struct axbVec_s *, void*)) op_desc.func;
  op(x, y, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecAXPY(struct axbVec_s *y, const struct axbScalar_s *alpha, const struct axbVec_s *x)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);

  axbOpDescriptor_t op_desc = y->opBackend->op_table[AXB_OP_VEC_AXPY];
  axbStatus_t (*op)(struct axbVec_s *, const struct axbScalar_s *, const struct axbVec_s *, void*) = (axbStatus_t (*)(struct axbVec_s *, const struct axbScalar_s *, const struct axbVec_s *, void*)) op_desc.func;
  op(y, alpha, x, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecAYPX(struct axbVec_s *y, const struct axbScalar_s *alpha, const struct axbVec_s *x)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);

  axbOpDescriptor_t op_desc = y->opBackend->op_table[AXB_OP_VEC_AYPX];
  axbStatus_t (*op)(struct axbVec_s *, const struct axbScalar_s *, const struct axbVec_s *, void*) = (axbStatus_t (*)(struct axbVec_s *, const struct axbScalar_s *, const struct axbVec_s *, void*)) op_desc.func;
  op(y, alpha, x, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecAXPBYPCZ(struct axbVec_s *z, const struct axbScalar_s *alpha, const struct axbScalar_s *beta, const struct axbScalar_s *gamma, const struct axbVec_s *x, const struct axbVec_s *y)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);
  AXB_VECTOR_INIT_CHECK(z);

  axbOpDescriptor_t op_desc = z->opBackend->op_table[AXB_OP_VEC_AXPBYPCZ];
  axbStatus_t (*op)(struct axbVec_s *, const struct axbScalar_s *, const struct axbScalar_s *, const struct axbScalar_s *, const struct axbVec_s *, const struct axbVec_s *, void*) = (axbStatus_t (*)(struct axbVec_s *, const struct axbScalar_s *, const struct axbScalar_s *, const struct axbScalar_s *, const struct axbVec_s *, const struct axbVec_s *, void*)) op_desc.func;
  op(z, alpha, beta, gamma, x, y, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecWAXPY(struct axbVec_s *w, const struct axbScalar_s *alpha, const struct axbVec_s *x, const struct axbVec_s *y)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);
  AXB_VECTOR_INIT_CHECK(w);

  axbOpDescriptor_t op_desc = w->opBackend->op_table[AXB_OP_VEC_WAXPY];
  axbStatus_t (*op)(struct axbVec_s *, const struct axbScalar_s *, const struct axbVec_s *, const struct axbVec_s *, void*) = (axbStatus_t (*)(struct axbVec_s *, const struct axbScalar_s *, const struct axbVec_s *, const struct axbVec_s *, void*)) op_desc.func;
  op(w, alpha, x, y, op_desc.func_data);
  return 0;
}

axbStatus_t axbVecMAXPY(struct axbVec_s *y, size_t num_vecs, const struct axbScalar_s * const *alpha, const struct axbVec_s * const *x)
{
  AXB_VECTOR_INIT_CHECK(y);
  for (size_t i=0; i<num_vecs; ++i) AXB_VECTOR_INIT_CHECK(x[i]);

  axbOpDescriptor_t op_desc = y->opBackend->op_table[AXB_OP_VEC_MAXPY];
  axbStatus_t (*op)(struct axbVec_s *, size_t, const struct axbScalar_s * const *, const struct axbVec_s * const *, void*) = (axbStatus_t (*)(struct axbVec_s *, size_t, const struct axbScalar_s * const *, const struct axbVec_s * const *, void*)) op_desc.func;
  op(y, num_vecs, alpha, x, op_desc.func_data);
  return 0;
}



axbStatus_t axbVecPointwiseMult(struct axbVec_s *w, const struct axbVec_s *x, const struct axbVec_s *y)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);
  AXB_VECTOR_INIT_CHECK(w);

  axbOpDescriptor_t op_desc = w->opBackend->op_table[AXB_OP_VEC_POINTWISEMULT];
  axbStatus_t (*op)(struct axbVec_s *, const struct axbVec_s *, const struct axbVec_s *, void*) = (axbStatus_t (*)(struct axbVec_s *, const struct axbVec_s *, const struct axbVec_s *, void*)) op_desc.func;
  op(w, x, y, op_desc.func_data);
  return 0;
}


axbStatus_t axbVecPointwiseDivide(struct axbVec_s *w, const struct axbVec_s *x, const struct axbVec_s *y)
{
  AXB_VECTOR_INIT_CHECK(x);
  AXB_VECTOR_INIT_CHECK(y);
  AXB_VECTOR_INIT_CHECK(w);

  axbOpDescriptor_t op_desc = w->opBackend->op_table[AXB_OP_VEC_POINTWISEDIV];
  axbStatus_t (*op)(struct axbVec_s *, const struct axbVec_s *, const struct axbVec_s *, void*) = (axbStatus_t (*)(struct axbVec_s *, const struct axbVec_s *, const struct axbVec_s *, void*)) op_desc.func;
  op(w, x, y, op_desc.func_data);
  return 0;
}


