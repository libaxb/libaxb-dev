#ifndef LIBAXB_BACKEND_OP_H_
#define LIBAXB_BACKEND_OP_H_

enum {
  // inplace operations:
  AXB_OP_VEC_SET = 0,
  AXB_OP_VEC_SQRTABS,
  AXB_OP_VEC_ZERO,
  AXB_OP_VEC_SCALE,

  // reduction operations
  AXB_OP_VEC_SUM,
  AXB_OP_VEC_DOT,
  AXB_OP_VEC_TDOT,
  AXB_OP_VEC_MDOT,
  AXB_OP_VEC_NORM1,
  AXB_OP_VEC_NORM2,
  AXB_OP_VEC_NORMINF,
  AXB_OP_VEC_DOTNORM2,
  AXB_OP_VEC_MAX,
  AXB_OP_VEC_MIN,

  // vector-vector operations
  AXB_OP_VEC_COPY,
  AXB_OP_VEC_SWAP,
  AXB_OP_VEC_AXPY,
  AXB_OP_VEC_AYPX,
  AXB_OP_VEC_AXPBYPCZ,
  AXB_OP_VEC_WAXPY,
  AXB_OP_VEC_MAXPY,
  AXB_OP_VEC_POINTWISEMULT,
  AXB_OP_VEC_POINTWISEDIV,

  //
  // MAT routines
  //
  AXB_OP_MAT_VEC,
  AXB_OP_MAT_TVEC,
  AXB_OP_MAT_MAT,
  AXB_OP_MAT_TRANS,

};


#endif
