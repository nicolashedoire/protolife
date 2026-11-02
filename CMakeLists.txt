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
  src/capsule.c
  src/capsule_shape.c
  src/diagnostic.c
  src/formula.c
  src/network.c
  src/ode.c
  src/protocol.c
  src/reactor.c
  src/rodas.c
  src/scene.c
  src/vesicles.c
  src/vesicles_balances.c
  src/vesicles_events.c)
add_library(protolife::protolife ALIAS protolife)
include(GNUInstallDirs)
target_include_directories(protolife PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(protolife PROPERTIES
  VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR} WINDOWS_EXPORT_ALL_SYMBOLS ON)
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

# Installing: the library, its headers, its data and a CMake package, so that another project can use
# find_package(protolife) and link protolife::protolife. A shared library carries SUNDIALS inside it, so it installs
# alone: configure with -DBUILD_SHARED_LIBS=ON to install. A static build installs the library and its headers only,
# since SUNDIALS's static libraries would have to be installed beside it.
include(CMakePackageConfigHelpers)
install(TARGETS protolife EXPORT protolifeTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/protolife DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.h")
install(FILES ${PROJECT_BINARY_DIR}/include/protolife/version.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/protolife)
install(DIRECTORY data/networks DESTINATION ${CMAKE_INSTALL_DATADIR}/protolife)
if(BUILD_SHARED_LIBS)
  set(PROTOLIFE_DATA_DIR ${CMAKE_INSTALL_DATADIR}/protolife)
  install(EXPORT protolifeTargets NAMESPACE protolife:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/protolife)
  configure_package_config_file(cmake/protolifeConfig.cmake.in ${PROJECT_BINARY_DIR}/protolifeConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/protolife PATH_VARS PROTOLIFE_DATA_DIR)
  write_basic_package_version_file(${PROJECT_BINARY_DIR}/protolifeConfigVersion.cmake COMPATIBILITY SameMinorVersion)
  install(FILES ${PROJECT_BINARY_DIR}/protolifeConfig.cmake ${PROJECT_BINARY_DIR}/protolifeConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/protolife)
endif()
