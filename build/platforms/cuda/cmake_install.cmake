# Install script for directory: /home/rc/vs-code-C/project2/openmm_image/platforms/cuda

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
  if(EXISTS "$ENV{DESTDIR}/home/rc/miniconda3/envs/openmm8.0/lib/plugins/libImagePluginCUDA.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/rc/miniconda3/envs/openmm8.0/lib/plugins/libImagePluginCUDA.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/rc/miniconda3/envs/openmm8.0/lib/plugins/libImagePluginCUDA.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/rc/miniconda3/envs/openmm8.0/lib/plugins/libImagePluginCUDA.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/rc/miniconda3/envs/openmm8.0/lib/plugins" TYPE SHARED_LIBRARY FILES "/home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda/libImagePluginCUDA.so")
  if(EXISTS "$ENV{DESTDIR}/home/rc/miniconda3/envs/openmm8.0/lib/plugins/libImagePluginCUDA.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/rc/miniconda3/envs/openmm8.0/lib/plugins/libImagePluginCUDA.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/rc/miniconda3/envs/openmm8.0/lib/plugins/libImagePluginCUDA.so"
         OLD_RPATH "/home/rc/miniconda3/envs/openmm8.0/lib:/home/rc/miniconda3/envs/openmm8.0/lib/plugins:/home/rc/vs-code-C/project2/openmm_image/build:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/rc/miniconda3/envs/openmm8.0/lib/plugins/libImagePluginCUDA.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda/tests/cmake_install.cmake")
endif()

