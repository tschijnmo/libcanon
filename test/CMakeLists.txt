#
# Build dependency Google test
#

set(BUILD_GTEST ON CACHE BOOL "Enable gtest") 
set(BUILD_GMOCK OFF CACHE BOOL "Disable gmock")
add_subdirectory(googletest EXCLUDE_FROM_ALL)

#
# Preparation for the actual tests
#

include(CTest)
function(add_simple_test name)
    add_executable("${name}" "${name}.cpp")
    target_include_directories("${name}" PRIVATE 
        "${libcanon_SOURCE_DIR}/include"
        "${gtest_SOURCE_DIR}/include")
    target_link_libraries("${name}" gtest gtest_main)
    add_test(NAME "${name}" COMMAND "${name}")
endfunction()

