#----------------------------------------------------------------
# Generated CMake target import file for configuration "MinSizeRel".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "SUNDIALS::core_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::core_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::core_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_core_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::core_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::core_static "${_IMPORT_PREFIX}/lib/sundials_core_static.lib" )

# Import target "SUNDIALS::core_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::core_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::core_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_core.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_core.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::core_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::core_shared "${_IMPORT_PREFIX}/lib/sundials_core.lib" "${_IMPORT_PREFIX}/bin/sundials_core.dll" )

# Import target "SUNDIALS::nvecserial_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::nvecserial_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::nvecserial_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_nvecserial_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::nvecserial_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::nvecserial_static "${_IMPORT_PREFIX}/lib/sundials_nvecserial_static.lib" )

# Import target "SUNDIALS::nvecserial_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::nvecserial_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::nvecserial_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_nvecserial.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_nvecserial.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::nvecserial_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::nvecserial_shared "${_IMPORT_PREFIX}/lib/sundials_nvecserial.lib" "${_IMPORT_PREFIX}/bin/sundials_nvecserial.dll" )

# Import target "SUNDIALS::nvecmanyvector_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::nvecmanyvector_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::nvecmanyvector_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_nvecmanyvector_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::nvecmanyvector_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::nvecmanyvector_static "${_IMPORT_PREFIX}/lib/sundials_nvecmanyvector_static.lib" )

# Import target "SUNDIALS::nvecmanyvector_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::nvecmanyvector_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::nvecmanyvector_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_nvecmanyvector.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_nvecmanyvector.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::nvecmanyvector_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::nvecmanyvector_shared "${_IMPORT_PREFIX}/lib/sundials_nvecmanyvector.lib" "${_IMPORT_PREFIX}/bin/sundials_nvecmanyvector.dll" )

# Import target "SUNDIALS::sunmatrixband_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunmatrixband_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunmatrixband_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunmatrixband_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunmatrixband_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunmatrixband_static "${_IMPORT_PREFIX}/lib/sundials_sunmatrixband_static.lib" )

# Import target "SUNDIALS::sunmatrixband_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunmatrixband_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunmatrixband_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunmatrixband.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_sunmatrixband.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunmatrixband_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunmatrixband_shared "${_IMPORT_PREFIX}/lib/sundials_sunmatrixband.lib" "${_IMPORT_PREFIX}/bin/sundials_sunmatrixband.dll" )

# Import target "SUNDIALS::sunmatrixdense_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunmatrixdense_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunmatrixdense_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunmatrixdense_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunmatrixdense_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunmatrixdense_static "${_IMPORT_PREFIX}/lib/sundials_sunmatrixdense_static.lib" )

# Import target "SUNDIALS::sunmatrixdense_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunmatrixdense_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunmatrixdense_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunmatrixdense.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_sunmatrixdense.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunmatrixdense_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunmatrixdense_shared "${_IMPORT_PREFIX}/lib/sundials_sunmatrixdense.lib" "${_IMPORT_PREFIX}/bin/sundials_sunmatrixdense.dll" )

# Import target "SUNDIALS::sunmatrixsparse_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunmatrixsparse_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunmatrixsparse_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunmatrixsparse_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunmatrixsparse_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunmatrixsparse_static "${_IMPORT_PREFIX}/lib/sundials_sunmatrixsparse_static.lib" )

# Import target "SUNDIALS::sunmatrixsparse_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunmatrixsparse_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunmatrixsparse_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunmatrixsparse.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_sunmatrixsparse.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunmatrixsparse_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunmatrixsparse_shared "${_IMPORT_PREFIX}/lib/sundials_sunmatrixsparse.lib" "${_IMPORT_PREFIX}/bin/sundials_sunmatrixsparse.dll" )

# Import target "SUNDIALS::sunlinsolband_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunlinsolband_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunlinsolband_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunlinsolband_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolband_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolband_static "${_IMPORT_PREFIX}/lib/sundials_sunlinsolband_static.lib" )

# Import target "SUNDIALS::sunlinsolband_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunlinsolband_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunlinsolband_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunlinsolband.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_sunlinsolband.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolband_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolband_shared "${_IMPORT_PREFIX}/lib/sundials_sunlinsolband.lib" "${_IMPORT_PREFIX}/bin/sundials_sunlinsolband.dll" )

# Import target "SUNDIALS::sunlinsoldense_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunlinsoldense_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunlinsoldense_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunlinsoldense_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsoldense_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsoldense_static "${_IMPORT_PREFIX}/lib/sundials_sunlinsoldense_static.lib" )

# Import target "SUNDIALS::sunlinsoldense_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunlinsoldense_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunlinsoldense_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunlinsoldense.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_sunlinsoldense.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsoldense_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsoldense_shared "${_IMPORT_PREFIX}/lib/sundials_sunlinsoldense.lib" "${_IMPORT_PREFIX}/bin/sundials_sunlinsoldense.dll" )

# Import target "SUNDIALS::sunlinsolpcg_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunlinsolpcg_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunlinsolpcg_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunlinsolpcg_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolpcg_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolpcg_static "${_IMPORT_PREFIX}/lib/sundials_sunlinsolpcg_static.lib" )

# Import target "SUNDIALS::sunlinsolpcg_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunlinsolpcg_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunlinsolpcg_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunlinsolpcg.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_sunlinsolpcg.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolpcg_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolpcg_shared "${_IMPORT_PREFIX}/lib/sundials_sunlinsolpcg.lib" "${_IMPORT_PREFIX}/bin/sundials_sunlinsolpcg.dll" )

# Import target "SUNDIALS::sunlinsolspbcgs_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunlinsolspbcgs_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunlinsolspbcgs_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunlinsolspbcgs_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolspbcgs_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolspbcgs_static "${_IMPORT_PREFIX}/lib/sundials_sunlinsolspbcgs_static.lib" )

# Import target "SUNDIALS::sunlinsolspbcgs_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunlinsolspbcgs_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunlinsolspbcgs_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunlinsolspbcgs.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_sunlinsolspbcgs.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolspbcgs_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolspbcgs_shared "${_IMPORT_PREFIX}/lib/sundials_sunlinsolspbcgs.lib" "${_IMPORT_PREFIX}/bin/sundials_sunlinsolspbcgs.dll" )

# Import target "SUNDIALS::sunlinsolspfgmr_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunlinsolspfgmr_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunlinsolspfgmr_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunlinsolspfgmr_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolspfgmr_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolspfgmr_static "${_IMPORT_PREFIX}/lib/sundials_sunlinsolspfgmr_static.lib" )

# Import target "SUNDIALS::sunlinsolspfgmr_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunlinsolspfgmr_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunlinsolspfgmr_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunlinsolspfgmr.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_sunlinsolspfgmr.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolspfgmr_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolspfgmr_shared "${_IMPORT_PREFIX}/lib/sundials_sunlinsolspfgmr.lib" "${_IMPORT_PREFIX}/bin/sundials_sunlinsolspfgmr.dll" )

# Import target "SUNDIALS::sunlinsolspgmr_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunlinsolspgmr_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunlinsolspgmr_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunlinsolspgmr_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolspgmr_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolspgmr_static "${_IMPORT_PREFIX}/lib/sundials_sunlinsolspgmr_static.lib" )

# Import target "SUNDIALS::sunlinsolspgmr_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunlinsolspgmr_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunlinsolspgmr_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunlinsolspgmr.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_sunlinsolspgmr.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolspgmr_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolspgmr_shared "${_IMPORT_PREFIX}/lib/sundials_sunlinsolspgmr.lib" "${_IMPORT_PREFIX}/bin/sundials_sunlinsolspgmr.dll" )

# Import target "SUNDIALS::sunlinsolsptfqmr_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunlinsolsptfqmr_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunlinsolsptfqmr_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunlinsolsptfqmr_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolsptfqmr_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolsptfqmr_static "${_IMPORT_PREFIX}/lib/sundials_sunlinsolsptfqmr_static.lib" )

# Import target "SUNDIALS::sunlinsolsptfqmr_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunlinsolsptfqmr_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunlinsolsptfqmr_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunlinsolsptfqmr.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_sunlinsolsptfqmr.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunlinsolsptfqmr_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunlinsolsptfqmr_shared "${_IMPORT_PREFIX}/lib/sundials_sunlinsolsptfqmr.lib" "${_IMPORT_PREFIX}/bin/sundials_sunlinsolsptfqmr.dll" )

# Import target "SUNDIALS::sunnonlinsolnewton_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunnonlinsolnewton_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunnonlinsolnewton_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunnonlinsolnewton_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunnonlinsolnewton_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunnonlinsolnewton_static "${_IMPORT_PREFIX}/lib/sundials_sunnonlinsolnewton_static.lib" )

# Import target "SUNDIALS::sunnonlinsolnewton_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunnonlinsolnewton_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunnonlinsolnewton_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunnonlinsolnewton.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_sunnonlinsolnewton.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunnonlinsolnewton_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunnonlinsolnewton_shared "${_IMPORT_PREFIX}/lib/sundials_sunnonlinsolnewton.lib" "${_IMPORT_PREFIX}/bin/sundials_sunnonlinsolnewton.dll" )

# Import target "SUNDIALS::sunnonlinsolfixedpoint_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunnonlinsolfixedpoint_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunnonlinsolfixedpoint_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunnonlinsolfixedpoint_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunnonlinsolfixedpoint_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunnonlinsolfixedpoint_static "${_IMPORT_PREFIX}/lib/sundials_sunnonlinsolfixedpoint_static.lib" )

# Import target "SUNDIALS::sunnonlinsolfixedpoint_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::sunnonlinsolfixedpoint_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::sunnonlinsolfixedpoint_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_sunnonlinsolfixedpoint.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_sunnonlinsolfixedpoint.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::sunnonlinsolfixedpoint_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::sunnonlinsolfixedpoint_shared "${_IMPORT_PREFIX}/lib/sundials_sunnonlinsolfixedpoint.lib" "${_IMPORT_PREFIX}/bin/sundials_sunnonlinsolfixedpoint.dll" )

# Import target "SUNDIALS::arkode_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::arkode_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::arkode_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_arkode_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::arkode_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::arkode_static "${_IMPORT_PREFIX}/lib/sundials_arkode_static.lib" )

# Import target "SUNDIALS::arkode_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::arkode_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::arkode_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_arkode.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_arkode.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::arkode_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::arkode_shared "${_IMPORT_PREFIX}/lib/sundials_arkode.lib" "${_IMPORT_PREFIX}/bin/sundials_arkode.dll" )

# Import target "SUNDIALS::cvode_static" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::cvode_static APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::cvode_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "C"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_cvode_static.lib"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::cvode_static )
list(APPEND _cmake_import_check_files_for_SUNDIALS::cvode_static "${_IMPORT_PREFIX}/lib/sundials_cvode_static.lib" )

# Import target "SUNDIALS::cvode_shared" for configuration "MinSizeRel"
set_property(TARGET SUNDIALS::cvode_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(SUNDIALS::cvode_shared PROPERTIES
  IMPORTED_IMPLIB_MINSIZEREL "${_IMPORT_PREFIX}/lib/sundials_cvode.lib"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/sundials_cvode.dll"
  )

list(APPEND _cmake_import_check_targets SUNDIALS::cvode_shared )
list(APPEND _cmake_import_check_files_for_SUNDIALS::cvode_shared "${_IMPORT_PREFIX}/lib/sundials_cvode.lib" "${_IMPORT_PREFIX}/bin/sundials_cvode.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
