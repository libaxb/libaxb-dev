
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)

if (ENABLE_CUDA)
  cuda_add_library(libaxb  backend/mem/host.c
                           backend/mem/cuda.cu
                           backend/op/host.c
                           backend/op/cuda.cu
                           backend/backend.c
                           vec.c
                           scalar.c
                           general.c)
ELSE(ENABLE_CUDA)
  add_library(libaxb  backend/mem/host.c
                      backend/mem/cuda.c
                      backend/op/host.c
                      backend/op/cuda.c
                      backend/backend.c
                      vec.c
                      scalar.c
                      general.c)
ENDIF(ENABLE_CUDA)
