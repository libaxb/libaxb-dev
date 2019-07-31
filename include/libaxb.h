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
 AXB_INT_32,
 AXB_INT_64,
 AXB_REAL_FLOAT,
 AXB_REAL_DOUBLE,
 AXB_COMPLEX_FLOAT,
 AXB_COMPLEX_DOUBLE
}                           axbDataType_t;
typedef int                 axbStatus_t;
struct axbHandle_s;
typedef struct axbHandle_s *axbHandle_t;


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
axbStatus_t axbInit(axbHandle_t *handle);


/** @brief Destroys the handle used by libaxb. Usually the last function to call from libaxb.
 *
 * @param handle    The handle obtained via axbInit() that should be destroyed.
 */
axbStatus_t axbFinalize(axbHandle_t handle);



///
////////// Error Handling
///

axbStatus_t axbErrorGetName(axbStatus_t error, char *name, int name_size);
axbStatus_t axbErrorGetString(axbStatus_t error, char *string, int string_size);



///
////////// Memory Backends
///
struct axbMemBackend_s;
typedef struct axbMemBackend_s    *axbMemBackend_t;

axbStatus_t axbMemBackendCreate(axbMemBackend_t *mem);
axbStatus_t axbMemBackendRegister(axbHandle_t handle, axbMemBackend_t mem);
axbStatus_t axbMemBackendSetName(axbMemBackend_t backend, const char *name);
axbStatus_t axbMemBackendGetName(axbMemBackend_t backend, const char **name);

/** @brief Returns all the available memory backends.
 *
 * @param handle     The libaxb context handle.
 * @param mem        Pointer to multiple MemBackEnd handles.
 * @param mem_size   In: Maximum number of elements to be written to `mem`. Out: Number of actual elements in `mem`.
*  @return          Returns a success- or error-code. @see axbErrorGetName(), axbErrorGetString()
 */
axbStatus_t axbMemBackendGetAll(axbHandle_t handle, axbMemBackend_t **mem, int *mem_size);
axbStatus_t axbMemBackendGetByName(axbHandle_t handle, axbMemBackend_t *mem, const char *name);

axbStatus_t axbMemBackendSetMalloc(axbMemBackend_t mem, axbStatus_t (*func)(void **, size_t, void *));
axbStatus_t axbMemBackendSetFree(axbMemBackend_t mem, axbStatus_t (*func)(void *, void *));

axbStatus_t axbMemBackendMalloc(axbMemBackend_t mem, size_t num_bytes, void **ptr);
axbStatus_t axbMemBackendFree(axbMemBackend_t mem, void *ptr);

axbStatus_t axbMemBackendSetCopyIn(axbMemBackend_t mem, axbStatus_t (*func)(void *, axbDataType_t, void *, axbDataType_t, size_t, void *));
axbStatus_t axbMemBackendSetCopyOut(axbMemBackend_t mem, axbStatus_t (*func)(void *, axbDataType_t, void *, axbDataType_t, size_t, void *));

axbStatus_t axbMemBackendCopyIn(axbMemBackend_t mem, void *src, axbDataType_t src_type, void *dest, axbDataType_t dest_type, size_t n);
axbStatus_t axbMemBackendCopyOut(axbMemBackend_t mem, void *src, axbDataType_t src_type, void *dest, axbDataType_t dest_type, size_t n);

axbStatus_t axbMemBackendSetDestroy(axbMemBackend_t mem, axbStatus_t (*func)(void*));
axbStatus_t axbMemBackendDestroy(axbMemBackend_t mem);


///
////////// Backends for operations
///
struct axbOpBackend_s;
typedef struct axbOpBackend_s     *axbOpBackend_t;
typedef int                        axbOperationID_t;

axbStatus_t axbOpBackendCreate(axbOpBackend_t *ops);
axbStatus_t axbOpBackendRegister(axbHandle_t handle, axbOpBackend_t ops);
axbStatus_t axbOpBackendSetName(axbOpBackend_t ops, const char *name);
axbStatus_t axbOpBackendGetName(axbOpBackend_t ops, const char **name);

/** @brief Adds a new worker routine.
 *
 * @param ops      The OpBackEnd handle for which the operation should be set.
 * @param op_name  The name of the operation. Maximum length is 32 characters.
 * @param op_func  Function pointer for the respective operation. The actual function signature is operation-dependent, so func is casted internally to the correct signature.
 * @param op_data  Pointer to auxiliary data that will be passed to func when called.
 * @param op_id    [Out] Identifier for the operation that can be used to quickly retrieve an operation later on.
 * @return         Returns a success- or error-code. @see axbErrorGetName(), axbErrorGetString()
 */
axbStatus_t axbOpBackendAddOperation(axbOpBackend_t ops, const char *op_name, axbStatus_t (*op_func)(void), void *op_data, axbOperationID_t *op_id);


/** @brief Returns all the available operations backends.
 *
 * @param handle     The libaxb context handle.
 * @param ops        Pointer to multiple operation backend handles.
 * @param ops_size   In: Maximum number of elements in ops, Out: Actual number of handles in `ops`.
*/
axbStatus_t axbOpBackendGetAll(axbHandle_t handle, axbOpBackend_t **ops, int *ops_size);
axbStatus_t axbOpBackendGetByName(axbHandle_t handle, axbOpBackend_t *ops, const char *name);

axbStatus_t axbOpBackendSetDestroy(axbOpBackend_t ops, axbStatus_t (*func)(void*));
axbStatus_t axbOpBackendDestroy(axbOpBackend_t ops);

///
////////// Scalar
///

/** @brief Opaque handle to a scalar */
struct axbScalar_s;
typedef struct axbScalar_s *axbScalar_t;

axbStatus_t axbScalarCreateBegin(axbHandle_t handle, axbScalar_t *scalar);
axbStatus_t axbScalarSetDataType(axbScalar_t scalar, axbDataType_t datatype);
axbStatus_t axbScalarSetBackend(axbScalar_t scalar, axbMemBackend_t mem);
axbStatus_t axbScalarCreateEnd(axbScalar_t scalar);

axbStatus_t axbScalarSetValue(axbScalar_t scalar, void *value, axbDataType_t value_datatype);

// convenience routine?
axbStatus_t axbScalarCreate(axbHandle_t handle, axbScalar_t *scalar, void *value, axbDataType_t datatype, axbMemBackend_t mem);

axbStatus_t axbScalarDestroy(axbScalar_t scalar);


///
////////// Vectors
///

/** @brief Opaque handle to a vector (dense or sparse) */
struct axbVec_s;
typedef struct axbVec_s *axbVec_t;

axbStatus_t axbVecCreateBegin(axbHandle_t handle, axbVec_t *vec);

axbStatus_t axbVecSetSize(axbVec_t vec, size_t size);
axbStatus_t axbVecGetSize(axbVec_t vec, size_t *size);

axbStatus_t axbVecSetDataType(axbVec_t vec, axbDataType_t datatype);
axbStatus_t axbVecGetDataType(axbVec_t vec, axbDataType_t *datatype);

axbStatus_t axbVecSetMemBackend(axbVec_t vec, axbMemBackend_t backend);
axbStatus_t axbVecGetMemBackend(axbVec_t vec, axbMemBackend_t *backend);

axbStatus_t axbVecSetOpBackend(axbVec_t vec, axbOpBackend_t backend);
axbStatus_t axbVecGetOpBackend(axbVec_t vec, axbOpBackend_t *backend);


axbStatus_t axbVecCreateEnd(axbVec_t vec);

axbStatus_t axbVecSetName(axbVec_t vec, const char *name);
axbStatus_t axbVecGetName(axbVec_t vec, const char **name);

axbStatus_t axbVecSetValues(axbVec_t vec, void *values, axbDataType_t values_datatype);
axbStatus_t axbVecGetValues(axbVec_t vec, void *values, axbDataType_t values_datatype);


// operations

/** @brief y = alpha * x + y */
axbStatus_t axbVecAXPY(axbVec_t y, axbScalar_t alpha, axbVec_t x);
// more to follow


axbStatus_t axbVecDestroy(axbVec_t vec);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
