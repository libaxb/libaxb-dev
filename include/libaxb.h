/**
 * libaxb - Library providing a single interface to many well-tuned GPU libraries
 *
 * (C) 2019, libaxb developers
 *
 * License: MIT license. See file LICENSE for details.
*/


typedef enum {
 AXB_INT_32,
 AXB_INT_64,
 AXB_REAL_FLOAT,
 AXB_REAL_DOUBLE,
 AXB_COMPLEX_FLOAT,
 AXB_COMPLEX_DOUBLE
}                          axbDataType_t;
typedef int                axbStatus_t;
typedef struct axbHandle_s axbHandle_t;


///////// General helper functions


/** @brief Entry point to libaxb. This is usually the first function from libaxb to be called. Creates an opaque handle where all helper data is stored.
*
*  @param handle    Pointer to the handle to initialize.
*  @return          Returns a success- or error-code. @see axbGetErrorName(), axbGetErrorString()
*/
axbStatus_t axbInit(axbHandle_t *handle);


/** @brief Destroys the handle used by libaxb. Usually the last function to call from libaxb.
 *
 * @param handle    The handle obtained via axbInit() that should be destroyed.
 */
axbStatus_t axbFinalize(axbHandle_t handle);



////////// Backends

typedef struct axbBackend_s    axbBackEnd_t;
typedef char *                 axbOperationID_t;

axbStatus_t axbBackendRegister(axbHandle_t handle, axbBackEnd_t *backend);
axbStatus_t axbBackendSetName(axbBackEnd_t backend, const char *name);
axbStatus_t axbBackendSetOperation(axbBackEnd_t backend, const axbOperationID_t op_identifier, void (*func)(void));



////////// Scalar

/** @brief Opaque handle to a scalar */
typedef struct axbScalar_s axbScalar_t;

axbStatus_t axbScalarCreateBegin(axbHandle_t handle, axbScalar_t *scalar);
axbStatus_t axbScalarSetDataType(axbScalar_t scalar, axbDataType_t datatype);
axbStatus_t axbScalarSetBackend(axbScalar_t scalar, axbBackEnd_t backend);
axbStatus_t axbScalarCreateEnd(axbScalar_t scalar);

axbStatus_t axbScalarSetValue(axbScalar_t scalar, void *value, axbDataType_t value_datatype);

// convenience routine?
axbStatus_t axbScalarCreate(axbHandle_t handle, axbScalar_t *scalar, void *value, axbDataType_t datatype, axbBackEnd_t backend);




////////// Vectors

/** @brief Opaque handle to a vector (dense or sparse) */
typedef struct axbVec_s axbVec_t;

axbStatus_t axbVecCreateBegin(axbVec_t *vec);

axbStatus_t axbVecSetSize(axbVec_t vec, int size);
axbStatus_t axbVecGetSize(axbVec_t vec, int *size);

axbStatus_t axbVecSetDataType(axbVec_t vec, axbDataType_t datatype);
axbStatus_t axbVecGetDataType(axbVec_t vec, axbDataType_t *datatype);

axbStatus_t axbVecSetBackend(axbVec_t vec, axbBackEnd_t backend);
axbStatus_t axbVecGetBackend(axbVec_t vec, axbBackEnd_t *backend);

axbStatus_t axbVecCreateEnd(axbVec_t vec);

axbStatus_t axbVecSetName(axbVec_t vec, const char *name);
axbStatus_t axbVecGetName(axbVec_t vec, const char **name);

axbStatus_t axbVecSetValues(axbVec_t vec, void *, axbDataType_t value_datatype);
axbStatus_t axbVecGetValues(axbVec_t vec, void *, axbDataType_t value_datatype);

axbStatus_t axbVecSetValues(axbVec_t vec, void *, axbDataType_t value_datatype);
axbStatus_t axbVecGetValues(axbVec_t vec, void *, axbDataType_t value_datatype);


// operations

/** @brief y = alpha * x + y */
axbStatus_t axbVecAXPY(axbVec_t y, axbScalar_t alpha, axbVec_t x);
// more to follow

