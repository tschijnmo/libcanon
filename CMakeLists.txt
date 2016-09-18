cmake_minimum_required(VERSION 3.6)
project(libcanon VERSION 0.1.0 LANGUAGES CXX)

# OPTIONS
include(GNUInstallDirs)
option(BUILD_TESTS "Build unit tests" ON)
option(BUILD_DOCS "Build the API documentation" ON)

# The main components.
add_subdirectory(include)

if (BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif ()

if (BUILD_DOCS)
    add_subdirectory(doc)
endif ()

# Generation of the CMake package files.
include(CMakePackageConfigHelpers)
set(CONFIG_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake")
set(CONFIG_FILE "${CMAKE_CURRENT_BINARY_DIR}/libcanon-config.cmake")
set(
        CONFIGVERSION_FILE
        "${CMAKE_CURRENT_BINARY_DIR}/libcanon-config-version.cmake"
)
configure_package_config_file("cmake/config.cmake.in" ${CONFIG_FILE}
        INSTALL_DESTINATION ${CONFIG_DESTINATION})
write_basic_package_version_file(${CONFIGVERSION_FILE}
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion)
install(EXPORT libcanon_exports FILE libcanon-targets.cmake
        NAMESPACE libcanon DESTINATION ${CONFIG_DESTINATION})
install(FILES ${CONFIG_FILE} ${CONFIGVERSION_FILE}
        DESTINATION ${CONFIG_DESTINATION})

# Old pkg-config support.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libcanon.pc.in
        ${CMAKE_CURRENT_BINARY_DIR}/libcanon.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libcanon.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)