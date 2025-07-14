#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "AMReX-Hydro::amrex_hydro_api" for configuration "Debug"
set_property(TARGET AMReX-Hydro::amrex_hydro_api APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(AMReX-Hydro::amrex_hydro_api PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/lib/amrex_hydro_api.lib"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/amrex_hydro_api.dll"
  )

list(APPEND _cmake_import_check_targets AMReX-Hydro::amrex_hydro_api )
list(APPEND _cmake_import_check_files_for_AMReX-Hydro::amrex_hydro_api "${_IMPORT_PREFIX}/lib/amrex_hydro_api.lib" "${_IMPORT_PREFIX}/bin/amrex_hydro_api.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
