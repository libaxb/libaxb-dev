
#include <string.h>
#include <stdio.h>
#include "libaxb.h"
#include "libaxb/general.h"
#include "libaxb/backend/op.h"


axbStatus_t axbMatCreateBegin(axbHandle_t handle, axbMat_t *mat) { (void)handle; (void)mat; return 0; }

axbStatus_t axbMatSetSizes(axbMat_t mat, size_t num_rows, size_t num_cols) { (void)mat; (void)num_rows; (void)num_cols; return 0; }
axbStatus_t axbMatGetSizes(axbMat_t mat, const size_t *num_rows, const size_t *num_cols) { (void)mat; (void)num_rows; (void)num_cols; return 0; }

axbStatus_t axbMatSetCSRIndexTypes(axbMat_t mat, axbDataType_t row_type, axbDataType_t col_type) { (void)mat; (void)row_type; (void)col_type; return 0; }
axbStatus_t axbMatGetCSRIndexTypes(axbMat_t mat, const axbDataType_t *row_type, const axbDataType_t *col_type) { (void)mat; (void)row_type; (void)col_type; return 0; }

axbStatus_t axbMatSetDataType(axbMat_t mat, axbDataType_t datatype) { (void)mat; (void)datatype; return 0; }
axbStatus_t axbMatGetDataType(axbMat_t mat, const axbDataType_t *datatype) { (void)mat; (void)datatype; return 0; }

axbStatus_t axbMatSetMemBackend(axbMat_t mat, axbMemBackend_t backend) { (void)mat; (void)backend; return 0; }
axbStatus_t axbMatGetMemBackend(axbMat_t mat, const axbMemBackend_t *backend) { (void)mat; (void)backend; return 0; }

axbStatus_t axbMatSetOpBackend(axbMat_t mat, axbOpBackend_t backend) { (void)mat; (void)backend; return 0; }
axbStatus_t axbMatGetOpBackend(axbMat_t mat, const axbOpBackend_t *backend) { (void)mat; (void)backend; return 0; }


axbStatus_t axbMatSetStorageType(axbMat_t mat, axbMatStorage_t storage) { (void)mat; (void)storage; return 0; }

axbStatus_t axbMatGetStorageType(axbMat_t mat, axbMatStorage_t *storage) { (void)mat; (void)storage; return 0; }

axbStatus_t axbMatCreateEnd(axbMat_t mat) { (void)mat; return 0; }

axbStatus_t axbMatSetName(axbMat_t mat, const char *name) { (void)mat; (void)name; return 0; }
axbStatus_t axbMatGetName(axbMat_t mat, const char *name, size_t max_name_size) { (void)mat; (void)name; (void)max_name_size; return 0; }

axbStatus_t axbMatSetValuesDense(axbMat_t mat, void *values, axbDataType_t values_datatype) { (void)mat; (void)values; (void)values_datatype; return 0; }

axbStatus_t axbMatGetValuesDense(axbMat_t mat, void *values, axbDataType_t values_datatype) { (void)mat; (void)values; (void)values_datatype; return 0; }

axbStatus_t axbMatGetNonzerosSize(axbMat_t mat, size_t *num_nonzeros) { (void)mat; (void)num_nonzeros; return 0; }

axbStatus_t axbMatSetValuesCSR(axbMat_t mat,
                               void *row_markers, axbDataType_t row_markers_datatype,
                               void *col_indices, axbDataType_t col_indices_datatype,
                               void *values, axbDataType_t values_datatype,
                               size_t num_values) { (void)mat; (void)row_markers; (void)row_markers_datatype; (void)col_indices; (void)col_indices_datatype; (void)values; (void)values_datatype; (void)num_values; return 0; }

axbStatus_t axbMatGetValuesCSR(axbMat_t mat,
                               void *row_markers, axbDataType_t row_markers_datatype,
                               void *col_indices, axbDataType_t col_indices_datatype,
                               void *values, axbDataType_t values_datatype) { (void)mat; (void)row_markers; (void)row_markers_datatype; (void)col_indices; (void)col_indices_datatype; (void)values; (void)values_datatype; return 0; }

axbStatus_t axbMatDestroy(axbMat_t mat) { (void)mat; return 0; }



// operations

axbStatus_t axbMatVec(axbMat_t A, axbVec_t x, axbVec_t Ax) { (void)A; (void)x; (void)Ax; return 0; }

axbStatus_t axbMatTVec(axbMat_t A, axbVec_t x, axbVec_t ATx) { (void)A; (void)x; (void)ATx; return 0; }

axbStatus_t axbMatMat(axbMat_t A, axbMat_t B, axbMat_t *AB) { (void)A; (void)B; (void)AB; return 0; }

axbStatus_t axbMatTrans(axbMat_t A, axbMat_t *AT) { (void)A; (void)AT; return 0; }

