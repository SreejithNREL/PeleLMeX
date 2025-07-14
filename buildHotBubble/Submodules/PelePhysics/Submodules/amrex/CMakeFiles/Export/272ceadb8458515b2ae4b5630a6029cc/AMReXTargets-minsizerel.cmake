#----------------------------------------------------------------
# Generated CMake target import file for configuration "MinSizeRel".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "AMReX::amrex_2d" for configuration "MinSizeRel"
set_property(TARGET AMReX::amrex_2d APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(AMReX::amrex_2d PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "CXX"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/amrex_2d.lib"
  )

list(APPEND _cmake_import_check_targets AMReX::amrex_2d )
list(APPEND _cmake_import_check_files_for_AMReX::amrex_2d "${_IMPORT_PREFIX}/lib/amrex_2d.lib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
