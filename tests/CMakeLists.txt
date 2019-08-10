include(CTest)
include_directories(${PROJECT_SOURCE_DIR}/include)

include_directories(${PROJECT_SOURCE_DIR}/tests/src)

ADD_EXECUTABLE(test-vec-ops src/vec-ops.c)
target_link_libraries(test-vec-ops   axb m)
ADD_TEST(vec-ops-host test-vec-ops "host")


ADD_EXECUTABLE(test-mat-sparse-ops src/mat-sparse-ops.c)
target_link_libraries(test-mat-sparse-ops   axb m)
ADD_TEST(mat-sparse-ops-host test-mat-sparse-ops "host")

if (ENABLE_CUDA)
  ADD_TEST(vec-ops-CUDA test-vec-ops "CUDA")
  ADD_TEST(mat-sparse-ops-CUDA test-vec-ops "CUDA")
endif(ENABLE_CUDA)

if (ENABLE_OPENCL)
  ADD_TEST(vec-ops-OpenCL test-vec-ops "OpenCL")
  ADD_TEST(mat-sparse-ops-OpenCL test-vec-ops "OpenCL")
endif(ENABLE_OPENCL)
