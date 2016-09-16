cmake_minimum_required(VERSION 3.6)
project(libcanon VERSION 0.1.0 LANGUAGES CXX)

# OPTIONS
include(GNUInstallDirs)
option(BUILD_TESTS "Build unit tests" ON)
option(BUILD_DOCS "Build the API documentation" ON)

include_directories(include/libcanon)

add_subdirectory(include)

if (BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif ()

if (BUILD_DOCS)
    add_subdirectory(doc)
endif ()
