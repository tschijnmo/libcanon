#
# Build dependency Google test
#

set(BUILD_GTEST ON CACHE BOOL "Enable gtest") 
set(BUILD_GMOCK OFF CACHE BOOL "Disable gmock")
add_subdirectory(googletest EXCLUDE_FROM_ALL)

