include(CTest)
include_directories(${PROJECT_SOURCE_DIR}/include)

ADD_EXECUTABLE(test-vec-axpy src/vec-axpy.c)
target_link_libraries(test-vec-axpy   libaxb)

ADD_TEST(vec-axpy-host test-vec-axpy "host")

if (ENABLE_CUDA)
  ADD_TEST(vec-axpy-CUDA test-vec-axpy "CUDA")
endif(ENABLE_CUDA)

if (ENABLE_OPENCL)
  ADD_TEST(vec-axpy-OpenCL test-vec-axpy "OpenCL")
endif(ENABLE_OPENCL)
