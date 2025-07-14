# Install script for directory: E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/src/sundials

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/PeleLMeX")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  MESSAGE("
Install shared components
")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "E:/PeleLMeX_Dev/buildHotBubble/bin/Debug/sundials_core_static.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "E:/PeleLMeX_Dev/buildHotBubble/bin/Release/sundials_core_static.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "E:/PeleLMeX_Dev/buildHotBubble/bin/MinSizeRel/sundials_core_static.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "E:/PeleLMeX_Dev/buildHotBubble/bin/RelWithDebInfo/sundials_core_static.lib")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY OPTIONAL FILES "E:/PeleLMeX_Dev/buildHotBubble/bin/Debug/sundials_core.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY OPTIONAL FILES "E:/PeleLMeX_Dev/buildHotBubble/bin/Release/sundials_core.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY OPTIONAL FILES "E:/PeleLMeX_Dev/buildHotBubble/bin/MinSizeRel/sundials_core.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY OPTIONAL FILES "E:/PeleLMeX_Dev/buildHotBubble/bin/RelWithDebInfo/sundials_core.lib")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE SHARED_LIBRARY FILES "E:/PeleLMeX_Dev/buildHotBubble/bin/Debug/sundials_core.dll")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE SHARED_LIBRARY FILES "E:/PeleLMeX_Dev/buildHotBubble/bin/Release/sundials_core.dll")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE SHARED_LIBRARY FILES "E:/PeleLMeX_Dev/buildHotBubble/bin/MinSizeRel/sundials_core.dll")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE SHARED_LIBRARY FILES "E:/PeleLMeX_Dev/buildHotBubble/bin/RelWithDebInfo/sundials_core.dll")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sundials" TYPE FILE FILES
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_adaptcontroller.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_band.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_base.hpp"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_context.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_context.hpp"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_convertibleto.hpp"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_core.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_core.hpp"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_dense.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_direct.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_errors.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_futils.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_iterative.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_linearsolver.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_linearsolver.hpp"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_logger.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_math.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_matrix.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_matrix.hpp"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_memory.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_memory.hpp"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_mpi_types.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_nonlinearsolver.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_nonlinearsolver.hpp"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_nvector.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_nvector.hpp"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_profiler.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_profiler.hpp"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_stepper.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_types_deprecated.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_types.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/sundials_version.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sundials/priv" TYPE FILE FILES
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/priv/sundials_context_impl.h"
    "E:/PeleLMeX_Dev/Submodules/PelePhysics/Submodules/sundials/include/sundials/priv/sundials_errors_impl.h"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "E:/PeleLMeX_Dev/buildHotBubble/Submodules/PelePhysics/Submodules/sundials/src/sundials/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
