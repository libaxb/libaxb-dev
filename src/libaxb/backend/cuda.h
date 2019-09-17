#ifndef LIBAXB_BACKEND_CUDA_H_
#define LIBAXB_BACKEND_CUDA_H_


#include "libaxb.h"

axbStatus_t axbMemBackendRegister_CUDA(struct axbHandle_s *handle);
axbStatus_t axbOpBackendRegister_CUDA(struct axbHandle_s *handle);

#endif
