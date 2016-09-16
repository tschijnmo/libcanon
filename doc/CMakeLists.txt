include(FindDoxygen)

if (NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen is required for building documentation")
endif ()

set(doxygen_file "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in" ${doxygen_file} @ONLY)

add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxygen_file}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Building the API documentation"
        VERBATIM
        )
