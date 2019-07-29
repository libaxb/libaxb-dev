#ifndef LIBAXB_BACKEND_CUDA_H_
#define LIBAXB_BACKEND_CUDA_H_


#include "libaxb.h"

axbStatus_t axbMemBackendRegister_CUDA(axbHandle_t handle);
axbStatus_t axbOpBackendRegister_CUDA(axbHandle_t handle);

#endif
