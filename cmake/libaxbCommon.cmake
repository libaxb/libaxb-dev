
include(CTest)
include(CMakeDependentOption)
include(AddCCompilerFlagIfSupported)
include(AddCLinkerFlagIfSupported)

# Installation directories
##########################

set(INSTALL_INCLUDE_DIR include CACHE PATH
   "Installation directory for headers")
if(WIN32 AND NOT CYGWIN)
   set(DEF_INSTALL_CMAKE_DIR CMake)
else()
   set(DEF_INSTALL_CMAKE_DIR lib/cmake/libaxb)
endif()
set(INSTALL_CMAKE_DIR "${DEF_INSTALL_CMAKE_DIR}" CACHE PATH
   "Installation directory for CMake files")

if(NOT IS_ABSOLUTE "${INSTALL_CMAKE_DIR}")
   set(INSTALL_CMAKE_DIR "${CMAKE_INSTALL_PREFIX}/${INSTALL_CMAKE_DIR}")
endif()
file(RELATIVE_PATH CONF_REL_INSTALL_PREFIX "${INSTALL_CMAKE_DIR}"
   "${CMAKE_INSTALL_PREFIX}")
if(NOT IS_ABSOLUTE "${INSTALL_INCLUDE_DIR}")
   set(INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/${INSTALL_INCLUDE_DIR}")
endif()
file(RELATIVE_PATH CONF_REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")

# User options
##############

option(ENABLE_CUDA "Use the CUDA backend" OFF)

option(BUILD_EXAMPLES "Build example programs" ON)

option(ENABLE_OPENCL "Use the OpenCL backend" OFF)

option(ENABLE_OPENMP "Use OpenMP acceleration" OFF)

option(ENABLE_ASAN "Build with address sanitizer if available" OFF)

option(ENABLE_PEDANTIC_FLAGS "Enable pedantic compiler flags (GCC and Clang only)" OFF)

mark_as_advanced(BOOSTPATH ENABLE_ASAN ENABLE_PEDANTIC_FLAGS)

# Find prerequisites
####################

if (ENABLE_CUDA)
   find_package(CUDA REQUIRED)
   set(CUDA_ARCH_FLAG "-arch=sm_50" CACHE STRING "Use one out of sm_30, sm_32, sm_35, sm_37, sm_50,  sm_52, sm_53, sm_60, ...")
   set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}" "${CUDA_ARCH_FLAG}" "-DLIBAXB_ENABLE_CUDA")
endif(ENABLE_CUDA)

if (ENABLE_OPENMP)
   find_package(OpenMP REQUIRED)
   set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}   -DLIBAXB_ENABLE_OPENMP")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DLIBAXB_ENABLE_OPENMP")
   set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS}    ${OpenMP_EXE_LINKER_FLAGS}")
   set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${OpenMP_MODULE_LINKER_FLAGS}")
   set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
   set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} ${OpenMP_STATIC_LINKER_FLAGS}")
endif(ENABLE_OPENMP)

if (ENABLE_ASAN)
  add_c_compiler_flag_if_supported("-fsanitize=address")
  add_c_linker_flag_if_supported("-fsanitize=address")
endif(ENABLE_ASAN)


if (ENABLE_OPENCL)
  find_package(OpenCL REQUIRED)
  set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}  ${OpenMP_C_FLAGS}    -DLIBAXB_ENABLE_OPENCL")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DLIBAXB_ENABLE_OPENCL")
endif(ENABLE_OPENCL)

# Set high warning level on GCC
if(ENABLE_PEDANTIC_FLAGS)
  set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -Wall -pedantic -Wextra -Wconversion")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -Wconversion")
endif()

# Export
########

configure_file(cmake/FindOpenCL.cmake
   "${PROJECT_BINARY_DIR}/FindOpenCL.cmake" COPYONLY)

configure_file(cmake/libaxbConfig.cmake.in
   "${PROJECT_BINARY_DIR}/libaxbConfig.cmake" @ONLY)

configure_file(cmake/libaxbConfigVersion.cmake.in
   "${PROJECT_BINARY_DIR}/libaxbConfigVersion.cmake" @ONLY)

if (CMAKE_MINOR_VERSION GREATER 6)  # export(PACKAGE ...) introduced with CMake 2.8.0
  export(PACKAGE libaxb)
endif()

# Library folders
#########

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Install
#########

install(FILES "${PROJECT_SOURCE_DIR}/include/libaxb.h"
   DESTINATION "${INSTALL_INCLUDE_DIR}" COMPONENT dev)

install(FILES
   "${PROJECT_BINARY_DIR}/FindOpenCL.cmake"
   "${PROJECT_BINARY_DIR}/libaxbConfig.cmake"
   "${PROJECT_BINARY_DIR}/libaxbConfigVersion.cmake"
   DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

# For out-of-the-box support on MacOS:
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  EXEC_PROGRAM(uname ARGS -v  OUTPUT_VARIABLE DARWIN_VERSION)
  STRING(REGEX MATCH "[0-9]+" DARWIN_VERSION ${DARWIN_VERSION})
  IF (DARWIN_VERSION GREATER 12)
    IF (ENABLE_CUDA)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")  # Mavericks and beyond need the new C++ STL with CUDA
    ENDIF()
  ENDIF()
  INCLUDE_DIRECTORIES("/opt/local/include")
  SET(CMAKE_EXE_LINKER_FLAGS "-framework OpenCL")
  set(CMAKE_MACOSX_RPATH 1) # Required for newer versions of CMake on MacOS X: http://www.kitware.com/blog/home/post/510
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
