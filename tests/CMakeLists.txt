include(CTest)
include_directories(${PROJECT_SOURCE_DIR}/include)

ADD_EXECUTABLE(test-vec-axpy src/vec-axpy.c)
target_link_libraries(test-vec-axpy   libaxb)
ADD_TEST(vec-axpy test-vec-axpy)
