include(CTest)
include_directories(${PROJECT_SOURCE_DIR}/include)

ADD_EXECUTABLE(test-vec-ops src/vec-ops.c)
target_link_libraries(test-vec-ops   libaxb m)

ADD_TEST(vec-ops-host test-vec-ops "host")

if (ENABLE_CUDA)
  ADD_TEST(vec-ops-CUDA test-vec-ops "CUDA")
endif(ENABLE_CUDA)

if (ENABLE_OPENCL)
  ADD_TEST(vec-ops-OpenCL test-vec-ops "OpenCL")
endif(ENABLE_OPENCL)
