
#include <string.h>
#include <stdio.h>
#include "libaxb.h"
#include "libaxb/general.h"
#include "libaxb/backend/op.h"


axbStatus_t axbMatCreateBegin(axbHandle_t handle, axbMat_t *mat)
{
  *mat = malloc(sizeof(struct axbMat_s));
  (*mat)->handle = handle;
  (*mat)->init = 119346;

  // set defaults:
  (*mat)->name = malloc(10);
  (*mat)->name[0] = 0;
  (*mat)->name_capacity = 10;

  (*mat)->rows = 0;
  (*mat)->cols = 0;
  (*mat)->nonzeros = 0;

  (*mat)->row_markers = NULL;
  (*mat)->row_markers_datatype = AXB_INT_32;

  (*mat)->col_indices = NULL;
  (*mat)->col_indices_datatype = AXB_INT_32;

  (*mat)->values = NULL;
  (*mat)->values_datatype = AXB_REAL_DOUBLE;

  (*mat)->memBackend = handle->memBackends[0];
  (*mat)->opBackend  = handle->opBackends[0];

  (*mat)->storage_type = AXB_STORAGE_CSR;

  return 0;
}

axbStatus_t axbMatSetSizes(axbMat_t mat, size_t num_rows, size_t num_cols)
{
  mat->rows = num_rows;
  mat->cols = num_cols;
  return 0;
}
axbStatus_t axbMatGetSizes(axbMat_t mat, size_t *num_rows, size_t *num_cols)
{
  *num_rows = mat->rows;
  *num_cols = mat->cols;
  return 0;
}

axbStatus_t axbMatSetCSRIndexTypes(axbMat_t mat, axbDataType_t row_type, axbDataType_t col_type)
{
  mat->row_markers_datatype = row_type;
  mat->col_indices_datatype = col_type;
  return 0;
}
axbStatus_t axbMatGetCSRIndexTypes(axbMat_t mat, axbDataType_t *row_type, axbDataType_t *col_type)
{
  *row_type = mat->row_markers_datatype;
  *col_type = mat->col_indices_datatype;
  return 0;
}

axbStatus_t axbMatSetDataType(axbMat_t mat, axbDataType_t datatype)
{
  mat->values_datatype = datatype;
  return 0;
}
axbStatus_t axbMatGetDataType(axbMat_t mat, axbDataType_t *datatype)
{
  *datatype = mat->values_datatype;
  return 0;
}

axbStatus_t axbMatSetMemBackend(axbMat_t mat, axbMemBackend_t backend)
{
  mat->memBackend = backend;
  return 0;
}
axbStatus_t axbMatGetMemBackend(axbMat_t mat, axbMemBackend_t *backend)
{
  *backend = mat->memBackend;
  return 0;
}

axbStatus_t axbMatSetOpBackend(axbMat_t mat, axbOpBackend_t backend)
{
  mat->opBackend = backend;
  return 0;
}
axbStatus_t axbMatGetOpBackend(axbMat_t mat, axbOpBackend_t *backend)
{
  *backend = mat->opBackend;
  return 0;
}


axbStatus_t axbMatSetStorageType(axbMat_t mat, axbMatStorage_t storage)
{
  mat->storage_type = storage;
  return 0;
}

axbStatus_t axbMatGetStorageType(axbMat_t mat, axbMatStorage_t *storage)
{
  *storage = mat->storage_type;
  return 0;
}

axbStatus_t axbMatCreateEnd(axbMat_t mat)
{
  mat->init = 119347;  // TODO: better abstraction for these magic numbers

  if (mat->storage_type == AXB_STORAGE_DENSE) {
    return axbMemBackendMalloc(mat->memBackend, sizeof(double) * (mat->rows * mat->cols), &(mat->values));
  }
  return 0;
}

axbStatus_t axbMatSetName(axbMat_t mat, const char *name)
{
  size_t len = strlen(name);
  if (len > mat->name_capacity) {
    free(mat->name);
    mat->name = malloc(len + 2);
  }
  for (size_t i=0; i<=len; ++i) mat->name[i] = name[i];
  return 0;
}
axbStatus_t axbMatGetName(axbMat_t mat, const char **name)
{
  *name = mat->name;
  return 0;
}

axbStatus_t axbMatSetValuesDense(axbMat_t mat, void *values, axbDataType_t values_datatype)
{
  if (mat->init != 119347) {
    fprintf(stderr, "ERROR in %s: Matrix not yet initialized. Call to axbMatCreateEnd() required.\n", __func__);
    return mat->init;
  }

  axbStatus_t status;
  if (mat->storage_type == AXB_STORAGE_CSR) {

    status = axbMemBackendFree(mat->memBackend, mat->row_markers); AXB_ERRCHK(status);
    status = axbMemBackendFree(mat->memBackend, mat->col_indices); AXB_ERRCHK(status);
    status = axbMemBackendFree(mat->memBackend, mat->values); AXB_ERRCHK(status);

    status = axbMemBackendMalloc(mat->memBackend, sizeof(int)    * (mat->rows+1),           &(mat->row_markers)); AXB_ERRCHK(status);
    status = axbMemBackendMalloc(mat->memBackend, sizeof(int)    * (mat->rows * mat->cols), &(mat->col_indices)); AXB_ERRCHK(status);
    status = axbMemBackendMalloc(mat->memBackend, sizeof(double) * (mat->rows * mat->cols), &(mat->values)); AXB_ERRCHK(status);

    // create temporary array for rows and columns:
    int *tmp_rows = malloc(sizeof(int)    * (mat->rows+1));
    int *tmp_cols = malloc(sizeof(int)    * (mat->rows * mat->cols));

    tmp_rows[0] = 0;
    int index = 0;
    for (int i=0; i<(int)mat->rows; ++i) {
      for (int j=0; j<(int)mat->cols; ++j) {
        tmp_cols[index++] = j;
      }
      tmp_rows[i+1] = index;
    }

    status = axbMemBackendCopyIn(mat->memBackend, values, values_datatype, mat->row_markers, mat->row_markers_datatype, mat->rows + 1); AXB_ERRCHK(status);
    status = axbMemBackendCopyIn(mat->memBackend, values, values_datatype, mat->col_indices, mat->col_indices_datatype, mat->rows * mat->cols); AXB_ERRCHK(status);
    status = axbMemBackendCopyIn(mat->memBackend, values, values_datatype, mat->values,      mat->values_datatype,      mat->rows * mat->cols); AXB_ERRCHK(status);

    free(tmp_rows);
    free(tmp_cols);
  }
  else if (mat->storage_type == AXB_STORAGE_DENSE)
  {
    status = axbMemBackendCopyIn(mat->memBackend, values, values_datatype, mat->values, mat->values_datatype, mat->rows * mat->cols); AXB_ERRCHK(status);
  } else {
    // TODO: implement
    return 1;
  }

  return 0;
}

axbStatus_t axbMatGetValuesDense(axbMat_t mat, void *values, axbDataType_t values_datatype)
{
  if (mat->storage_type == AXB_STORAGE_DENSE) {
    return axbMemBackendCopyOut(mat->memBackend, mat->values, mat->values_datatype, values, values_datatype, mat->rows * mat->cols);
  } else {
    // TODO: implement
    return 1;
  }
  return 0;
}

axbStatus_t axbMatGetNonzerosSize(axbMat_t mat, size_t *num_nonzeros)
{
  if (mat->storage_type == AXB_STORAGE_DENSE) { // count the number of actual nonzeros
    size_t nnz = 0;
    for (size_t i=0; i<mat->rows; ++i)
      for (size_t j=0; j<mat->cols; ++j) {
        double val = ((double*)mat->values)[i*mat->cols + j];
        if (val < 0 || val > 0) ++nnz;
      }
    *num_nonzeros = nnz;
  } else {
    *num_nonzeros = mat->nonzeros;
  }
  return 0;
}

axbStatus_t axbMatSetValuesCSR(axbMat_t mat,
                               void *row_markers, axbDataType_t row_markers_datatype,
                               void *col_indices, axbDataType_t col_indices_datatype,
                               void *values, axbDataType_t values_datatype,
                               size_t num_values)
{
  axbStatus_t status;
  if (mat->storage_type == AXB_STORAGE_CSR) // just copy over
  {
    if (mat->nonzeros != num_values) {
      mat->nonzeros = num_values;
      status = axbMemBackendMalloc(mat->memBackend, sizeof(int)    * (mat->rows+1),   &(mat->row_markers)); AXB_ERRCHK(status);
      status = axbMemBackendMalloc(mat->memBackend, sizeof(int)    * (mat->nonzeros), &(mat->col_indices)); AXB_ERRCHK(status);
      status = axbMemBackendMalloc(mat->memBackend, sizeof(double) * (mat->nonzeros), &(mat->values)); AXB_ERRCHK(status);
    }
    status = axbMemBackendCopyIn(mat->memBackend, row_markers, row_markers_datatype, mat->row_markers, mat->row_markers_datatype, mat->rows + 1); AXB_ERRCHK(status);
    status = axbMemBackendCopyIn(mat->memBackend, col_indices, col_indices_datatype, mat->col_indices, mat->col_indices_datatype, num_values); AXB_ERRCHK(status);
    status = axbMemBackendCopyIn(mat->memBackend, values,      values_datatype,      mat->values,      mat->values_datatype,      num_values); AXB_ERRCHK(status);
  } else if (mat->storage_type == AXB_STORAGE_DENSE) {
    int *r = (int*)row_markers;
    int *c = (int*)col_indices;
    double *v = (double*)values;

    double *tmp_values = malloc(sizeof(double) * mat->rows * mat->cols);

    for (int i=0; i<(int)mat->rows; ++i) {
      for (int j=r[i]; j<r[i+1]; ++j) {
        tmp_values[(size_t)i * mat->cols + (size_t)c[j]] = v[j];
      }
    }

    status = axbMemBackendCopyIn(mat->memBackend, tmp_values, values_datatype, mat->values, values_datatype, mat->rows * mat->cols); AXB_ERRCHK(status);
  } else {
    assert( 1==2 && "axbMatSetValuesCSR() for compressed CSR not yet implemented!");
  }
  return 0;
}

axbStatus_t axbMatGetValuesCSR(axbMat_t mat,
                               void *row_markers, axbDataType_t row_markers_datatype,
                               void *col_indices, axbDataType_t col_indices_datatype,
                               void *values, axbDataType_t values_datatype)
{
  axbStatus_t status;
  if (mat->storage_type == AXB_STORAGE_DENSE) {
    int *r = row_markers;
    int *c = col_indices;
    double *v = values;

    double *values = mat->values;

    int index = 0;
    r[0] = 0;
    for (size_t i=0; i<mat->rows; ++i) {
      for (size_t j=0; j<mat->cols; ++j) {
        double val = values[i];
        if (val < 0 || val > 0) {
          c[index] = (int)j;
          v[index] = val;
          ++index;
        }
      }
      r[i+1] = index; // terminate row
    }
  } else if (mat->storage_type == AXB_STORAGE_CSR) {
    status = axbMemBackendCopyOut(mat->memBackend, mat->row_markers, mat->row_markers_datatype, row_markers, row_markers_datatype, mat->rows + 1); AXB_ERRCHK(status);
    status = axbMemBackendCopyOut(mat->memBackend, mat->col_indices, mat->col_indices_datatype, col_indices, col_indices_datatype, mat->nonzeros); AXB_ERRCHK(status);
    status = axbMemBackendCopyOut(mat->memBackend,      mat->values,      mat->values_datatype, values,      values_datatype,      mat->nonzeros); AXB_ERRCHK(status);
  } else {
    assert( 1==2 && "axbMatGetValuesCSR() for compressed CSR not yet implemented!");
  }
  return 0;
}

axbStatus_t axbMatDestroy(axbMat_t mat)
{
  // TODO: Check proper value of mat->init
  mat->init += 1;

  axbStatus_t status;
  status = axbMemBackendFree(mat->memBackend, mat->row_markers); AXB_ERRCHK(status);
  status = axbMemBackendFree(mat->memBackend, mat->col_indices); AXB_ERRCHK(status);
  status = axbMemBackendFree(mat->memBackend, mat->values); AXB_ERRCHK(status);

  free(mat->name);
  free(mat);
  return 0;
}



// operations

axbStatus_t axbMatVec(axbMat_t A, axbVec_t x, axbVec_t Ax)
{
  axbOpDescriptor_t op_desc = A->opBackend->op_table[AXB_OP_MAT_VEC];
  axbStatus_t (*op)(axbMat_t, axbVec_t, axbVec_t, void*) = (axbStatus_t (*)(axbMat_t, axbVec_t, axbVec_t, void*)) op_desc.func;
  op(A, x, Ax, op_desc.func_data);
  return 0;
}

axbStatus_t axbMatTVec(axbMat_t A, axbVec_t x, axbVec_t ATx)
{
  axbOpDescriptor_t op_desc = A->opBackend->op_table[AXB_OP_MAT_TVEC];
  axbStatus_t (*op)(axbMat_t, axbVec_t, axbVec_t, void*) = (axbStatus_t (*)(axbMat_t, axbVec_t, axbVec_t, void*)) op_desc.func;
  op(A, x, ATx, op_desc.func_data);
  return 0;
}

axbStatus_t axbMatMat(axbMat_t A, axbMat_t B, axbMat_t *AB)
{
  axbOpDescriptor_t op_desc = A->opBackend->op_table[AXB_OP_MAT_MAT];
  axbStatus_t (*op)(axbMat_t, axbMat_t, axbMat_t*, void*) = (axbStatus_t (*)(axbMat_t, axbMat_t, axbMat_t*, void*)) op_desc.func;
  op(A, B, AB, op_desc.func_data);
  return 0;
}

axbStatus_t axbMatTrans(axbMat_t A, axbMat_t *AT)
{
  axbOpDescriptor_t op_desc = A->opBackend->op_table[AXB_OP_MAT_TRANS];
  axbStatus_t (*op)(axbMat_t, axbMat_t*, void*) = (axbStatus_t (*)(axbMat_t, axbMat_t*, void*)) op_desc.func;
  op(A, AT, op_desc.func_data);
  return 0;
}

