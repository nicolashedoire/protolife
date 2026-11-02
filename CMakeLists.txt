cmake_minimum_required(VERSION 3.24)
project(protolife VERSION 0.0.1 DESCRIPTION "Building blocks for synthetic life, in simulation" LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(PROTOLIFE_TESTS "Build the tests" ${PROJECT_IS_TOP_LEVEL})
option(PROTOLIFE_SANITIZE "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)

include(cmake/sundials.cmake)

configure_file(include/protolife/version.h.in include/protolife/version.h @ONLY)

add_library(protolife
  src/diagnostic.c
  src/formula.c
  src/network.c
  src/ode.c
  src/reactor.c)
add_library(protolife::protolife ALIAS protolife)
target_include_directories(protolife PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(protolife PRIVATE SUNDIALS::cvode)
if(NOT MSVC)
  target_link_libraries(protolife PUBLIC m)
endif()

# Warnings for protolife's own code only, never for the dependencies.
add_library(protolife_warnings INTERFACE)
if(MSVC)
  target_compile_options(protolife_warnings INTERFACE /W4 /WX)
  target_compile_definitions(protolife_warnings INTERFACE _CRT_SECURE_NO_WARNINGS)
else()
  target_compile_options(protolife_warnings INTERFACE
    -Wall -Wextra -Wpedantic -Wshadow -Wstrict-prototypes -Werror)
endif()
target_link_libraries(protolife PRIVATE $<BUILD_INTERFACE:protolife_warnings>)

if(PROTOLIFE_SANITIZE AND NOT MSVC)
  target_compile_options(protolife PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(protolife PUBLIC -fsanitize=address,undefined)
endif()

if(PROTOLIFE_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

option(PROTOLIFE_EXAMPLES "Build the examples" ${PROJECT_IS_TOP_LEVEL})
if(PROTOLIFE_EXAMPLES)
  add_subdirectory(examples)
endif()
