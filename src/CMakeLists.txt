
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)

add_library(libaxpy  backend/mem/host.c
                     backend/op/host.c
                     backend/backend.c
                     vec.c
                     scalar.c
                     general.c)
