#ifndef LIBAXP_LIBAXB_H
#define LIBAXP_LIBAXB_H

/**
 * libaxb - Library providing a single interface to many well-tuned GPU libraries
 *
 * (C) 2019, libaxb developers
 *
 * License: MIT license. See file LICENSE for details.
*/

#include <stdlib.h>

typedef enum {
 AXB_INT_16,
 AXB_INT_32,
 AXB_INT_64,
 AXB_REAL_HALF,
 AXB_REAL_FLOAT,
 AXB_REAL_DOUBLE,
 AXB_COMPLEX_FLOAT,
 AXB_COMPLEX_DOUBLE
}                           axbDataType_t;
typedef int                 axbStatus_t;
struct axbHandle_s;


#ifdef __cplusplus
extern "C"
{
#endif

///
///////// General helper functions
///

/** @brief Entry point to libaxb. This is usually the first function from libaxb to be called. Creates an opaque handle where all helper data is stored.
*
*  @param handle    Pointer to the handle to initialize.
*  @return          Returns a success- or error-code. @see axbErrorGetName(), axbErrorGetString()
*/
axbStatus_t axbInit(struct axbHandle_s **handle);


/** @brief Destroys the handle used by libaxb. Usually the last function to call from libaxb.
 *
 * @param handle    The handle obtained via axbInit() that should be destroyed.
 */
axbStatus_t axbFinalize(struct axbHandle_s *handle);



///
////////// Error Handling
///

axbStatus_t axbErrorGetName(axbStatus_t error, char *name, int name_size);
axbStatus_t axbErrorGetString(axbStatus_t error, char *string, int string_size);



///
////////// Memory Backends
///
struct axbMemBackend_s;

axbStatus_t axbMemBackendCreate(struct axbMemBackend_s **mem);
axbStatus_t axbMemBackendRegister(struct axbHandle_s *handle, struct axbMemBackend_s *mem);
axbStatus_t axbMemBackendSetName(struct axbMemBackend_s *backend, const char *name);
axbStatus_t axbMemBackendGetName(struct axbMemBackend_s *backend, const char **name);

/** @brief Returns all the available memory backends.
 *
 * @param handle     The libaxb context handle.
 * @param mem        Pointer to multiple MemBackEnd handles.
 * @param mem_size   In: Maximum number of elements to be written to `mem`. Out: Number of actual elements in `mem`.
*  @return          Returns a success- or error-code. @see axbErrorGetName(), axbErrorGetString()
 */
axbStatus_t axbMemBackendGetAll(struct axbHandle_s *handle, struct axbMemBackend_s ***mem, size_t *mem_size);
axbStatus_t axbMemBackendGetByName(struct axbHandle_s *handle, struct axbMemBackend_s **mem, const char *name);

axbStatus_t axbMemBackendSetMalloc(struct axbMemBackend_s *mem, axbStatus_t (*func)(void **, size_t, void *));
axbStatus_t axbMemBackendSetFree(struct axbMemBackend_s *mem, axbStatus_t (*func)(void *, void *));

axbStatus_t axbMemBackendMalloc(struct axbMemBackend_s *mem, size_t num_bytes, void **ptr);
axbStatus_t axbMemBackendFree(struct axbMemBackend_s *mem, void *ptr);

axbStatus_t axbMemBackendSetCopyIn(struct axbMemBackend_s *mem, axbStatus_t (*func)(void *, axbDataType_t, void *, axbDataType_t, size_t, void *));
axbStatus_t axbMemBackendSetCopyOut(struct axbMemBackend_s *mem, axbStatus_t (*func)(void *, axbDataType_t, void *, axbDataType_t, size_t, void *));

axbStatus_t axbMemBackendCopyIn(struct axbMemBackend_s *mem, void *src, axbDataType_t src_type, void *dest, axbDataType_t dest_type, size_t n);
axbStatus_t axbMemBackendCopyOut(struct axbMemBackend_s *mem, void *src, axbDataType_t src_type, void *dest, axbDataType_t dest_type, size_t n);

axbStatus_t axbMemBackendSetDestroy(struct axbMemBackend_s *mem, axbStatus_t (*func)(void*));
axbStatus_t axbMemBackendDestroy(struct axbMemBackend_s *mem);


///
////////// Backends for operations
///
struct axbOpBackend_s;
typedef int axbOperationID_t;

axbStatus_t axbOpBackendCreate(struct axbOpBackend_s **ops);
axbStatus_t axbOpBackendRegister(struct axbHandle_s *handle, struct axbOpBackend_s *ops);
axbStatus_t axbOpBackendSetName(struct axbOpBackend_s *ops, const char *name);
axbStatus_t axbOpBackendGetName(struct axbOpBackend_s *ops, const char **name);

/** @brief Adds a new worker routine.
 *
 * @param ops      The OpBackEnd handle for which the operation should be set.
 * @param op_name  The name of the operation. Maximum length is 32 characters.
 * @param op_func  Function pointer for the respective operation. The actual function signature is operation-dependent, so func is casted internally to the correct signature.
 * @param op_data  Pointer to auxiliary data that will be passed to func when called.
 * @param op_id    [Out] Identifier for the operation that can be used to quickly retrieve an operation later on.
 * @return         Returns a success- or error-code. @see axbErrorGetName(), axbErrorGetString()
 */
axbStatus_t axbOpBackendAddOperation(struct axbOpBackend_s *ops, const char *op_name, axbStatus_t (*op_func)(void), void *op_data, axbOperationID_t *op_id);


/** @brief Returns all the available operations backends.
 *
 * @param handle     The libaxb context handle.
 * @param ops        Pointer to multiple operation backend handles.
 * @param ops_size   In: Maximum number of elements in ops, Out: Actual number of handles in `ops`.
*/
axbStatus_t axbOpBackendGetAll(struct axbHandle_s *handle, struct axbOpBackend_s ***ops, size_t *ops_size);
axbStatus_t axbOpBackendGetByName(struct axbHandle_s *handle, struct axbOpBackend_s **ops, const char *name);

axbStatus_t axbOpBackendSetDestroy(struct axbOpBackend_s *ops, axbStatus_t (*func)(void*));
axbStatus_t axbOpBackendDestroy(struct axbOpBackend_s *ops);

/*
 * Scalar
 */

/** @brief Opaque handle to a scalar */
struct axbScalar_s;

axbStatus_t axbScalarCreateBegin(struct axbHandle_s *handle, struct axbScalar_s **scalar);
axbStatus_t axbScalarSetDataType(struct axbScalar_s *scalar, axbDataType_t datatype);
axbStatus_t axbScalarSetBackend(struct axbScalar_s *scalar, struct axbMemBackend_s *mem);
axbStatus_t axbScalarCreateEnd(struct axbScalar_s *scalar);

axbStatus_t axbScalarSetValue(struct axbScalar_s *scalar, void *value, axbDataType_t value_datatype);
axbStatus_t axbScalarGetValue(const struct axbScalar_s *scalar, void *value, axbDataType_t value_datatype);

// convenience routine?
axbStatus_t axbScalarCreate(struct axbHandle_s *handle, struct axbScalar_s **scalar, void *value, axbDataType_t datatype, struct axbMemBackend_s *mem);

axbStatus_t axbScalarDestroy(struct axbScalar_s *scalar);


/*
 * Vectors
 */

/** @brief Opaque handle to a vector (dense or sparse) */
struct axbVec_s;

axbStatus_t axbVecCreateBegin(struct axbHandle_s * handle, struct axbVec_s **vec);

axbStatus_t axbVecSetSize(struct axbVec_s *vec, size_t size);
axbStatus_t axbVecGetSize(const struct axbVec_s *vec, size_t *size);

axbStatus_t axbVecSetDataType(struct axbVec_s *vec, axbDataType_t datatype);
axbStatus_t axbVecGetDataType(const struct axbVec_s *vec, axbDataType_t *datatype);

axbStatus_t axbVecSetMemBackend(struct axbVec_s *vec, struct axbMemBackend_s *backend);
axbStatus_t axbVecGetMemBackend(const struct axbVec_s *vec, struct axbMemBackend_s **backend);

axbStatus_t axbVecSetOpBackend(struct axbVec_s *vec, struct axbOpBackend_s *backend);
axbStatus_t axbVecGetOpBackend(const struct axbVec_s *vec, struct axbOpBackend_s **backend);


axbStatus_t axbVecCreateEnd(struct axbVec_s *vec);

axbStatus_t axbVecSetName(struct axbVec_s *vec, const char *name);
axbStatus_t axbVecGetName(const struct axbVec_s *vec, const char **name);

axbStatus_t axbVecSetValues(struct axbVec_s *vec, void *values, axbDataType_t values_datatype);
axbStatus_t axbVecGetValues(const struct axbVec_s *vec, void *values, axbDataType_t values_datatype);

axbStatus_t axbVecDestroy(struct axbVec_s *vec);


// operations

// initialize vector:
axbStatus_t axbVecSet(struct axbVec_s *x, const struct axbScalar_s *value);
//TODO: axbStatus_t axbVecSetValues(struct axbVec_s *x, void *indices, size_t num_indices, axbDataType_t indices_datatype, void *values, size_t num_values, axbDataType_t values_datatype);

// in-place operations
/** @brief Replaces all entries of the vector with the sqrt of the absolute value */
axbStatus_t axbVecSqrtAbs(struct axbVec_s *x);
/** @brief Sets all vector entries to zero */
axbStatus_t axbVecZero(struct axbVec_s *x);
/** @brief Scales the vector x by a scalar factor alpha, i.e. x[i] <- alpha * x[i]  */
axbStatus_t axbVecScale(struct axbVec_s *x, const struct axbScalar_s *alpha);

// reduction operations:
/** @brief Computes the sum of entries in x. Usually faster than computing the inner product (x,1) with '1' being a vector of ones. */
axbStatus_t axbVecSum(const struct axbVec_s *x, struct axbScalar_s *sum);
/** @brief Computes the dot-product (x,y) = y^H x, where H denotes the conjugate transpose of y. */
axbStatus_t axbVecDot(const struct axbVec_s *x, const struct axbVec_s *y, struct axbScalar_s *dot);
/** @brief Computes the indefinite dot-product (x,y), i.e. no complex conjugation on either x or y */
axbStatus_t axbVecTDot(const struct axbVec_s *x, const struct axbVec_s *y, struct axbScalar_s *tdot);
/** @brief Computes the inner products (x, y[0]), (x, y[1]), ..., (x, y[num_vecs-1]) */
axbStatus_t axbVecMDot(const struct axbVec_s *x, size_t num_vecs, const struct axbVec_s **y, struct axbScalar_s **mdot);
/** @brief Computes the 1-norm of the vector x */
axbStatus_t axbVecNorm1(const struct axbVec_s *x, struct axbScalar_s *norm);
/** @brief Computes the 2-norm of the vector x */
axbStatus_t axbVecNorm2(const struct axbVec_s *x, struct axbScalar_s *norm);
/** @brief Computes the inf-norm of the vector x */
axbStatus_t axbVecNormInf(const struct axbVec_s *x, struct axbScalar_s *norm);
/** @brief Computes the inner product s*t and the 2-norm of t */
axbStatus_t axbVecDotNorm2(const struct axbVec_s *s, const struct axbVec_s *t, struct axbScalar_s *dot_st, struct axbScalar_s *norm_t);

/** @brief Computes the element with the maximum real part and its location
*
*   @param idx    The smallest index of an element with the maximum value
*/
axbStatus_t axbVecMax(const struct axbVec_s *x, size_t *idx, struct axbScalar_s *m);
/** @brief Computes the element with the minimum real part and its location
*
*   @param idx    The smallest index of an element with the minimum value
*/
axbStatus_t axbVecMin(const struct axbVec_s *x, size_t *idx, struct axbScalar_s *m);

// vector-vector operations:
/** @brief Assigns x to y */
axbStatus_t axbVecCopy(const struct axbVec_s *x, struct axbVec_s *y);
/** @brief Swaps the vectors x and y */
axbStatus_t axbVecSwap(struct axbVec_s *x, struct axbVec_s *y);

/** @brief y = alpha * x + y */
axbStatus_t axbVecAXPY(struct axbVec_s *y, const struct axbScalar_s *alpha, const struct axbVec_s *x);
/** @brief y = x + alpha * y */
axbStatus_t axbVecAYPX(struct axbVec_s *y, const struct axbScalar_s *alpha, const struct axbVec_s *x);
/** @brief z = alpha * x + beta * y + gamma * z */
axbStatus_t axbVecAXPBYPCZ(struct axbVec_s *z, const struct axbScalar_s *alpha, const struct axbScalar_s *beta, const struct axbScalar_s *gamma, const struct axbVec_s *x, const struct axbVec_s *y);
/** @brief w = alpha * x + y */
axbStatus_t axbVecWAXPY(struct axbVec_s *w, const struct axbScalar_s *alpha, const struct axbVec_s *x, const struct axbVec_s *y);
/** @brief y = y + \sum_i alpha[i] * x[i] */
axbStatus_t axbVecMAXPY(struct axbVec_s *y, size_t num_vecs, const struct axbScalar_s * const *alpha, const struct axbVec_s * const *x);

/** @brief w = x .* y  (component-wise multiplication) */
axbStatus_t axbVecPointwiseMult(struct axbVec_s *w, const struct axbVec_s *x, const struct axbVec_s *y);
/** @brief w = x ./ y  (component-wise division) */
axbStatus_t axbVecPointwiseDivide(struct axbVec_s *w, const struct axbVec_s *x, const struct axbVec_s *y);






/*
 * Matrices
 */

/** @brief Opaque handle to a matrix (dense or sparse) */
struct axbMat_s;

/** @brief Creates the memory representing the matrix object, starting initialization.
 *
 * Customize the matrix via calls to axbMatSetSizes(), axbMatSetCSRIndexTypes(), axbMatSetDataType(), axbMatSetMemBackend(), axbMatSetOpBackend(), axbMatSetStorageType()
 * A call to axbMatCreateEnd() will complete the initialization.
 *
 * @param handle     The surrounding axb environment
 * @param mat        Pointer to the matrix object to be created
 */
axbStatus_t axbMatCreateBegin(struct axbHandle_s *handle, struct axbMat_s **mat);

/** @brief Sets the matrix dimensions of the matrix
 *
 * @param mat      The matrix for which the size should be set
 * @param num_rows The number of rows of the matrix. Must be larger than zero.
 * @param num_cols The number of columns of the matrix. Must be larger than zero.
 */
axbStatus_t axbMatSetSizes(struct axbMat_s *mat, size_t num_rows, size_t num_cols);

/** @brief Returns the number of rows and columns of the matrix */
axbStatus_t axbMatGetSizes(const struct axbMat_s *mat, size_t *num_rows, size_t *num_cols);

/** @brief Specifies the index types to be used when the matrix is stored in a CSR format (i.e. AXB_STORAGE_CSR or AXB_STORAGE_COMPRESSED_CSR).
 *
 * Usually AXB_INT_32 is fine for the indices (up to ~2 billion nonzeros).
 * Larger problems should try with row_type of AXB_INT_64 first (while still benefitting from faster execution due to 32-bit column indices).
 * Only for matrices with a few rows and MANY columns it may be required to set col_type to AXB_INT_64.
 *
 * Must be called after axbMatCreateBegin() and before axbMatCreateEnd() for the particular matrix.
 *
 * @param mat         The matrix for which the index types should be set
 * @param row_type    The index type to use for the row-markers in the CSR format. Default: AXB_INT_32
 * @param col_type    The index type to use for the column indices in the CSR format. Default: AXB_INT_32
 */
axbStatus_t axbMatSetCSRIndexTypes(struct axbMat_s *mat, axbDataType_t row_type, axbDataType_t col_type);

/** @brief Returns the index data types used for the row-markers (row_type) and the column indices (col_type) for the CSR format. */
axbStatus_t axbMatGetCSRIndexTypes(const struct axbMat_s *mat, axbDataType_t *row_type, axbDataType_t *col_type);

/** @brief Sets the data type used for numerical entries in 'mat'. Default: AXB_REAL_DOUBLE
*
*  Must be called after axbMatCreateBegin() and before axbMatCreateEnd() for the particular matrix.
*/
axbStatus_t axbMatSetDataType(struct axbMat_s *mat, axbDataType_t datatype);

/** @brief Returns the data type used for the numerical entries in 'mat' */
axbStatus_t axbMatGetDataType(const struct axbMat_s *mat, axbDataType_t *datatype);

/** @brief Sets the memory backend to be used with the matrix.
 *
 *  Must be called after axbMatCreateBegin() and before axbMatCreateEnd() for the particular matrix.
 */
axbStatus_t axbMatSetMemBackend(struct axbMat_s *mat, struct axbMemBackend_s *backend);
/** @brief Returns the memory backend currently used by the matrix */
axbStatus_t axbMatGetMemBackend(const struct axbMat_s *mat, struct axbMemBackend_s **backend);

/** @brief Sets the operation backend used for running operations on the matrix */
axbStatus_t axbMatSetOpBackend(struct axbMat_s *mat, struct axbOpBackend_s *backend);

/** @brief Returns the operation backend currently in use for the matrix */
axbStatus_t axbMatGetOpBackend(const struct axbMat_s *mat, struct axbOpBackend_s **backend);

typedef enum {
 AXB_STORAGE_CSR = 0,           /// standard CSR format
 AXB_STORAGE_COMPRESSED_CSR,    /// compressed CSR format (suitable for matrices where most rows are zero)
 AXB_STORAGE_DENSE              /// row-major dense format
} axbMatStorage_t;

/** @brief Sets the internal storage format for the matrix.
 *
 * Must be called after axbMatCreateBegin() and before axbMatCreateEnd() for the particular matrix
 */
axbStatus_t axbMatSetStorageType(struct axbMat_s *mat, axbMatStorage_t storage);

/** @brief Returns the internal storage format used by the matrix */
axbStatus_t axbMatGetStorageType(const struct axbMat_s *mat, axbMatStorage_t *storage);

axbStatus_t axbMatCreateEnd(struct axbMat_s *mat);

/** @brief Sets the name for the matrix */
axbStatus_t axbMatSetName(struct axbMat_s *mat, const char *name);
/** @brief Returns the pointer to the name of the matrix */
axbStatus_t axbMatGetName(const struct axbMat_s *mat, const char **name);

/** @brief Populates the matrix with the provided values (row-major).
 *
 * @param mat               The matrix to be filled
 * @param values            Pointer to the values. Must be of size [number of rows] * [number of columns].
 * @param values_datatype   Data type descriptor for the values array (e.g. AXB_REAL_DOUBLE)
 */
axbStatus_t axbMatSetValuesDense(struct axbMat_s *mat, void *values, axbDataType_t values_datatype);

/** @brief Returns all values (including zeros) for the matrix in a row-major fashion.
 *
 * @param mat               The matrix for which the values should be returned
 * @param values            Pointer to the values. Must be at least of size [number of rows] * [number of columns].
 * @param values_datatype   Data type descriptor for the values array (e.g. AXB_REAL_DOUBLE)
 */
axbStatus_t axbMatGetValuesDense(const struct axbMat_s *mat, void *values, axbDataType_t values_datatype);

/** @brief Returns the number of nonzero entries in 'mat'
 *
 * @param mat           Matrix (either sparse or dense)
 * @param num_nonzeros  The number of nonzeros
 */
axbStatus_t axbMatGetNonzerosSize(const struct axbMat_s *mat, size_t *num_nonzeros);

/** @brief Sets the values of the matrix based on data in CSR format (sparse matrix)
 *
 *  Example for CSR for the matrix
 *    ( 1    0    2 )
 *    ( 3    0    0 )
 *    ( 0    4    5 )
 *  row_markers:   0,    2, 3,    5   [size num_rows + 1]
 *  col_indices:   0, 2, 0, 1, 2      [size num_nonzeros]
 *  values:        1, 2, 3, 4, 5      [size num_nonzeros]
 *
 * Can be used for dense matrices as well (but will not be efficient if there are no zeros in the matrix)
 *
 * @param row_markers    Entry points for each row into 'col_indices' and 'values'
 * @param col_indices    Column indices in ascending order (!)
 * @param values         The values to be set at the particular row and column
 */
axbStatus_t axbMatSetValuesCSR(struct axbMat_s *mat, void *row_markers, axbDataType_t row_markers_datatype, void *col_indices, axbDataType_t col_indices_datatype, void *values, axbDataType_t values_datatype, size_t num_values);

/** @brief Returns the values of the matrix in CSR format.
 *
 * Can be used for dense matrices as well (but will not be efficient if there are no zeros in the matrix)
 *
 * @see axbMatSetValuesCSR
 */
axbStatus_t axbMatGetValuesCSR(const struct axbMat_s *mat, void *row_markers, axbDataType_t row_markers_datatype, void *col_indices, axbDataType_t col_indices_datatype, void *values, axbDataType_t values_datatype);

/** @brief Destroys the matrix. */
axbStatus_t axbMatDestroy(struct axbMat_s *mat);



// operations

/** @brief Computes the matrix-vector product A * x for a matrix A and a vector x.
 *
 * @param A  The matrix
 * @param x  The vector to be multiplied with A
 * @param Ax The result vector. The size of Ax needs to be the same as the number of rows of A.
 */
axbStatus_t axbMatVec(const struct axbMat_s *A, const struct axbVec_s *x, struct axbVec_s *Ax);

/** @brief Computes the transposed matrix-vector product A^T * x for a matrix A and a vector x.
 *
 * @param A   The matrix
 * @param x   The vector to be multiplied with A
 * @param ATx The result vector. The size of ATx needs to be the same as the number of columns of A.
 */
axbStatus_t axbMatTVec(const struct axbMat_s *A, const struct axbVec_s *x, struct axbVec_s *ATx);

/** @brief Computes the matrix-matrix product A * B.
 *
 * @param A   The left factor. The number of columns in A must match the number of rows in B.
 * @param B   The right factor. The number of rows in B must match the number of columns in A.
 * @param AB  The result matrix. Will be newly created based on A and B.
 */
axbStatus_t axbMatMat(const struct axbMat_s *A, const struct axbMat_s *B, struct axbMat_s **AB);

/** @brief Computes the transpose A^T of A
 *
 * @param A   The matrix to be transposed
 * @param AT  The result of transposing A. Will be newly created based on A.
 */
axbStatus_t axbMatTrans(const struct axbMat_s *A, struct axbMat_s **AT);

// TODO: Add more operations



#ifdef __cplusplus
} // extern "C"
#endif

#endif
