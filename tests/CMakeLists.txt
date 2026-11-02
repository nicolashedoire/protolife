foreach(name ode network reactor vesicles rodas scene capsule)
  add_executable(test_${name} test_${name}.c)
  target_link_libraries(test_${name} PRIVATE protolife protolife_warnings)
  target_compile_definitions(test_${name} PRIVATE PROTOLIFE_DATA="${PROJECT_SOURCE_DIR}/data"
                                                  PROTOLIFE_TEST_DIR="${CMAKE_CURRENT_BINARY_DIR}")
  add_test(NAME ${name} COMMAND test_${name})
endforeach()
