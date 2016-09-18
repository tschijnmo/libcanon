# The phony licanon library for export.
add_library(headers INTERFACE)
target_include_directories(headers INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libcanon>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/libcanon>)
install(TARGETS headers EXPORT libcanon_exports)

# Installation of the real thing.
install(DIRECTORY libcanon
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        COMPONENT Devel)
