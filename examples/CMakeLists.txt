
include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(vec-axpy   vec-axpy.c)

target_link_libraries(vec-axpy  libaxpy)
