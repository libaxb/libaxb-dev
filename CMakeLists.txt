# Project setup
###############

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


project(libaxb)

set(VERSION_MAJOR 0)
set(VERSION_MINOR 1)
set(VERSION_PATCH 0)
set(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Add visibility of headers
# Necessary for Qt-Creator usage.
################
#file( GLOB_RECURSE MAKE_HEADERS_VISIBLE_SRC *.hpp *.h)
#add_custom_target( MAKE_HEADERS_VISIBLE SOURCES ${MAKE_HEADERS_VISIBLE_SRC} )


# Common options
################

include(libaxbCommon)

# Set release build:
IF(DEFINED CMAKE_BUILD_TYPE)
 SET (CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
ELSE()
 SET (CMAKE_BUILD_TYPE Release)
ENDIF()

# User options
##############

option(BUILD_DOXYGEN_DOCS "Build the Doxygen-generated API docs" OFF)

cmake_dependent_option(ENABLE_DIST "Enable the the distribution tarball target" ON "BUILD_DOXYGEN_DOCS" OFF)

mark_as_advanced(BUILD_DOXYGEN_DOCS ENABLE_DIST)

# Find prerequisites
####################

if(BUILD_DOXYGEN_DOCS)
   find_package(Doxygen REQUIRED)
endif()


# Compile options
#################

include_directories(BEFORE "${PROJECT_SOURCE_DIR}")


# Subdirectories
################

if(BUILD_EXAMPLES)
   add_subdirectory(examples)
endif()

if(BUILD_TESTING)
   INCLUDE(CTest)
   add_subdirectory(tests)
endif()

add_subdirectory(doc)

add_subdirectory(src)


