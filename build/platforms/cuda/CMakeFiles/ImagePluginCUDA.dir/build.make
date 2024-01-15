# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/rc/miniconda3/envs/openmm8.0/bin/cmake

# The command to remove a file.
RM = /home/rc/miniconda3/envs/openmm8.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/rc/vs-code-C/project2/openmm_image

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/rc/vs-code-C/project2/openmm_image/build

# Include any dependencies generated for this target.
include platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/compiler_depend.make

# Include the progress variables for this target.
include platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/progress.make

# Include the compile flags for this target's objects.
include platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/flags.make

platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/src/CudaImageKernelFactory.cpp.o: platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/flags.make
platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/src/CudaImageKernelFactory.cpp.o: /home/rc/vs-code-C/project2/openmm_image/platforms/cuda/src/CudaImageKernelFactory.cpp
platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/src/CudaImageKernelFactory.cpp.o: platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rc/vs-code-C/project2/openmm_image/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/src/CudaImageKernelFactory.cpp.o"
	cd /home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/src/CudaImageKernelFactory.cpp.o -MF CMakeFiles/ImagePluginCUDA.dir/src/CudaImageKernelFactory.cpp.o.d -o CMakeFiles/ImagePluginCUDA.dir/src/CudaImageKernelFactory.cpp.o -c /home/rc/vs-code-C/project2/openmm_image/platforms/cuda/src/CudaImageKernelFactory.cpp

platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/src/CudaImageKernelFactory.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ImagePluginCUDA.dir/src/CudaImageKernelFactory.cpp.i"
	cd /home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rc/vs-code-C/project2/openmm_image/platforms/cuda/src/CudaImageKernelFactory.cpp > CMakeFiles/ImagePluginCUDA.dir/src/CudaImageKernelFactory.cpp.i

platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/src/CudaImageKernelFactory.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ImagePluginCUDA.dir/src/CudaImageKernelFactory.cpp.s"
	cd /home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rc/vs-code-C/project2/openmm_image/platforms/cuda/src/CudaImageKernelFactory.cpp -o CMakeFiles/ImagePluginCUDA.dir/src/CudaImageKernelFactory.cpp.s

platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernels.cpp.o: platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/flags.make
platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernels.cpp.o: /home/rc/vs-code-C/project2/openmm_image/platforms/common/src/CommonImageKernels.cpp
platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernels.cpp.o: platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rc/vs-code-C/project2/openmm_image/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernels.cpp.o"
	cd /home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernels.cpp.o -MF CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernels.cpp.o.d -o CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernels.cpp.o -c /home/rc/vs-code-C/project2/openmm_image/platforms/common/src/CommonImageKernels.cpp

platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernels.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernels.cpp.i"
	cd /home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rc/vs-code-C/project2/openmm_image/platforms/common/src/CommonImageKernels.cpp > CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernels.cpp.i

platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernels.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernels.cpp.s"
	cd /home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rc/vs-code-C/project2/openmm_image/platforms/common/src/CommonImageKernels.cpp -o CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernels.cpp.s

platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernelSources.cpp.o: platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/flags.make
platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernelSources.cpp.o: platforms/common/src/CommonImageKernelSources.cpp
platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernelSources.cpp.o: platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rc/vs-code-C/project2/openmm_image/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernelSources.cpp.o"
	cd /home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernelSources.cpp.o -MF CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernelSources.cpp.o.d -o CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernelSources.cpp.o -c /home/rc/vs-code-C/project2/openmm_image/build/platforms/common/src/CommonImageKernelSources.cpp

platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernelSources.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernelSources.cpp.i"
	cd /home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rc/vs-code-C/project2/openmm_image/build/platforms/common/src/CommonImageKernelSources.cpp > CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernelSources.cpp.i

platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernelSources.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernelSources.cpp.s"
	cd /home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rc/vs-code-C/project2/openmm_image/build/platforms/common/src/CommonImageKernelSources.cpp -o CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernelSources.cpp.s

# Object files for target ImagePluginCUDA
ImagePluginCUDA_OBJECTS = \
"CMakeFiles/ImagePluginCUDA.dir/src/CudaImageKernelFactory.cpp.o" \
"CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernels.cpp.o" \
"CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernelSources.cpp.o"

# External object files for target ImagePluginCUDA
ImagePluginCUDA_EXTERNAL_OBJECTS =

platforms/cuda/libImagePluginCUDA.so: platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/src/CudaImageKernelFactory.cpp.o
platforms/cuda/libImagePluginCUDA.so: platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernels.cpp.o
platforms/cuda/libImagePluginCUDA.so: platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/__/common/src/CommonImageKernelSources.cpp.o
platforms/cuda/libImagePluginCUDA.so: platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/build.make
platforms/cuda/libImagePluginCUDA.so: /usr/local/cuda-10.2/lib64/libcudart_static.a
platforms/cuda/libImagePluginCUDA.so: /usr/lib/x86_64-linux-gnu/librt.so
platforms/cuda/libImagePluginCUDA.so: libOpenMMImage.so
platforms/cuda/libImagePluginCUDA.so: platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/rc/vs-code-C/project2/openmm_image/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX shared library libImagePluginCUDA.so"
	cd /home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ImagePluginCUDA.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/build: platforms/cuda/libImagePluginCUDA.so
.PHONY : platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/build

platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/clean:
	cd /home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda && $(CMAKE_COMMAND) -P CMakeFiles/ImagePluginCUDA.dir/cmake_clean.cmake
.PHONY : platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/clean

platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/depend:
	cd /home/rc/vs-code-C/project2/openmm_image/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/rc/vs-code-C/project2/openmm_image /home/rc/vs-code-C/project2/openmm_image/platforms/cuda /home/rc/vs-code-C/project2/openmm_image/build /home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda /home/rc/vs-code-C/project2/openmm_image/build/platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : platforms/cuda/CMakeFiles/ImagePluginCUDA.dir/depend

