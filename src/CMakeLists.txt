
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)

if (ENABLE_CUDA)
  cuda_add_library(axb     backend/mem/host.c
                           backend/mem/cuda.cu
                           backend/mem/opencl.c
                           backend/op/host.c
                           backend/op/cuda.cu
                           backend/op/opencl.c
                           backend/backend.c
                           mat.c
                           vec.c
                           scalar.c
                           general.c)
ELSE(ENABLE_CUDA)
  add_library(axb     backend/mem/host.c
                      backend/mem/cuda.c
                      backend/mem/opencl.c
                      backend/op/host.c
                      backend/op/cuda.c
                      backend/op/opencl.c
                      backend/backend.c
                      mat.c
                      vec.c
                      scalar.c
                      general.c)
ENDIF(ENABLE_CUDA)

IF (ENABLE_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIRS})
  target_link_libraries(axb ${OPENCL_LIBRARIES})
ENDIF(ENABLE_OPENCL)

target_link_libraries(axb m)

install(TARGETS axb
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
