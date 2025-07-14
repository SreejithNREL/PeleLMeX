#add_subdirectory(Plasma)
#add_subdirectory(Production)
#add_subdirectory(RegTests)
if(DEFINED PELELMEX_CASE)
    add_subdirectory(RegTests/${PELELMEX_CASE})
else()
    add_subdirectory(RegTests)  # or build all by default
endif()
#add_subdirectory(UnitTests)
