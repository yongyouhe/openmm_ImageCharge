# Install script for directory: /home/rc/vs-code-C/project2/openmm_image

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/rc/miniconda3/envs/openmm8.0")
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

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/rc/miniconda3/envs/openmm8.0/lib/libOpenMMImage.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/rc/miniconda3/envs/openmm8.0/lib/libOpenMMImage.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/rc/miniconda3/envs/openmm8.0/lib/libOpenMMImage.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/rc/miniconda3/envs/openmm8.0/lib/libOpenMMImage.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/rc/miniconda3/envs/openmm8.0/lib" TYPE SHARED_LIBRARY FILES "/home/rc/vs-code-C/project2/openmm_image/build/libOpenMMImage.so")
  if(EXISTS "$ENV{DESTDIR}/home/rc/miniconda3/envs/openmm8.0/lib/libOpenMMImage.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/rc/miniconda3/envs/openmm8.0/lib/libOpenMMImage.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/rc/miniconda3/envs/openmm8.0/lib/libOpenMMImage.so"
         OLD_RPATH "/home/rc/miniconda3/envs/openmm8.0/lib:/home/rc/miniconda3/envs/openmm8.0/lib/plugins:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/rc/miniconda3/envs/openmm8.0/lib/libOpenMMImage.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES "/home/rc/vs-code-C/project2/openmm_image/openmmapi/include/OpenMMImage.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/openmm" TYPE FILE FILES
    "/home/rc/vs-code-C/project2/openmm_image/openmmapi/include/openmm/ImageCustomIntegrator.h"
    "/home/rc/vs-code-C/project2/openmm_image/openmmapi/include/openmm/ImageIntegrator.h"
    "/home/rc/vs-code-C/project2/openmm_image/openmmapi/include/openmm/ImageKernels.h"
    "/home/rc/vs-code-C/project2/openmm_image/openmmapi/include/openmm/ImageLangevinIntegrator.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/openmm/internal" TYPE FILE FILES "/home/rc/vs-code-C/project2/openmm_image/openmmapi/include/openmm/internal/ImageCustomIntegratorUtilities.h")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/rc/vs-code-C/project2/openmm_image/build/serialization/tests/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/rc/vs-code-C/project2/openmm_image/build/platforms/common/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/rc/vs-code-C/project2/openmm_image/build/python/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/rc/vs-code-C/project2/openmm_image/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
