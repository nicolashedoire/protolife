foreach(name ode network reactor vesicles rodas scene capsule)
  add_executable(test_${name} test_${name}.c)
  target_link_libraries(test_${name} PRIVATE protolife protolife_warnings)
  target_compile_definitions(test_${name} PRIVATE PROTOLIFE_DATA="${PROJECT_SOURCE_DIR}/data"
                                                  PROTOLIFE_TEST_DIR="${CMAKE_CURRENT_BINARY_DIR}")
  add_test(NAME ${name} COMMAND test_${name})
endforeach()

# Step 8 shortened (test_selection.c): a cycle of the reproduction example's world in each of its three arms, about
# a minute each, labelled long so that slow builds (the sanitizers) can leave them out.
add_executable(test_selection test_selection.c ${PROJECT_SOURCE_DIR}/examples/reproduction_world.c)
target_include_directories(test_selection PRIVATE ${PROJECT_SOURCE_DIR}/examples)
target_link_libraries(test_selection PRIVATE protolife protolife_warnings)
target_compile_definitions(test_selection PRIVATE PROTOLIFE_DATA="${PROJECT_SOURCE_DIR}/data")
foreach(arm living twins stopped)
  add_test(NAME selection_${arm} COMMAND test_selection ${arm})
  set_tests_properties(selection_${arm} PROPERTIES LABELS long)
endforeach()
