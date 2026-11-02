# A project that uses an installed protolife, as anyone's would: built in CI to check the installation.
cmake_minimum_required(VERSION 3.24)
project(downstream LANGUAGES C)
find_package(protolife 0.0 REQUIRED)
add_executable(downstream main.c)
target_link_libraries(downstream PRIVATE protolife::protolife)
target_compile_definitions(downstream PRIVATE NETWORKS="${protolife_DATA_DIR}/networks")
