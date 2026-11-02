# SUNDIALS (BSD-3-Clause), built from its release with only CVODE, as a static library inside protolife.
# To build without network access, pass -DFETCHCONTENT_SOURCE_DIR_SUNDIALS=/path/to/sundials-7.9.0.
include(FetchContent)

set(SUNDIALS_ENABLE_CVODE ON CACHE BOOL "" FORCE)
foreach(package ARKODE CVODES IDA IDAS KINSOL)
  set(SUNDIALS_ENABLE_${package} OFF CACHE BOOL "" FORCE)
endforeach()
set(SUNDIALS_ENABLE_C_EXAMPLES OFF CACHE BOOL "" FORCE)
set(SUNDIALS_ENABLE_CXX_EXAMPLES OFF CACHE BOOL "" FORCE)
set(SUNDIALS_ENABLE_EXAMPLES_INSTALL OFF CACHE BOOL "" FORCE)
set(SUNDIALS_ENABLE_BENCHMARKS OFF CACHE BOOL "" FORCE)
# Failures reach protolife through the error handler and become a pl_diagnostic; nothing is printed.
set(SUNDIALS_LOGGING_LEVEL 0 CACHE STRING "" FORCE)
set(BUILD_STATIC_LIBS ON CACHE BOOL "" FORCE)
set(_protolife_shared ${BUILD_SHARED_LIBS})
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

FetchContent_Declare(sundials
  URL https://github.com/llnl/sundials/releases/download/v7.9.0/sundials-7.9.0.tar.gz
  URL_HASH SHA256=13f898a27b48fe3449483f9e438a800ed545abf93bc2e2ceec2d1e00ae8db5ef
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  EXCLUDE_FROM_ALL)
FetchContent_MakeAvailable(sundials)

set(BUILD_SHARED_LIBS ${_protolife_shared} CACHE BOOL "Build protolife as a shared library" FORCE)
