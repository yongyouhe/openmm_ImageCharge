#---------------------------------------------------
# OpenMM Image Plugin CUDA Platform
#----------------------------------------------------

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_SOURCE_SUBDIRS . ../common) # a list

# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(IMAGE_CUDA_LIBRARY_NAME ImagePluginCUDA)

SET(SHARED_TARGET ${IMAGE_CUDA_LIBRARY_NAME})

# 需要openmm的头文件.在openmm源码里的插件不需要加这个,因为是和openmm源码一起编译的,在40-48行添加了源码头文件
INCLUDE_DIRECTORIES(BEFORE "${OPENMM_DIR}/include/openmm/common" "${OPENMM_DIR}/include/openmm/reference" "${OPENMM_DIR}/include/openmm/cuda")

# These are all the places to search for header files which are
# to be part of the API.
SET(API_INCLUDE_DIRS)
FOREACH(subdir ${OPENMM_SOURCE_SUBDIRS})
    #append
    SET(API_INCLUDE_DIRS ${API_INCLUDE_DIRS}
                         ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include
                         ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include/internal)
ENDFOREACH(subdir)

# Locate header files.
# We'll need both *relative* path names, starting with their API_INCLUDE_DIRS,
# and absolute pathnames.
SET(API_REL_INCLUDE_FILES)   # start these out empty
SET(API_ABS_INCLUDE_FILES)

FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)	# returns full pathnames
    SET(API_ABS_INCLUDE_FILES ${API_ABS_INCLUDE_FILES} ${fullpaths})

    FOREACH(pathname ${fullpaths})
        GET_FILENAME_COMPONENT(filename ${pathname} NAME)
        SET(API_REL_INCLUDE_FILES ${API_REL_INCLUDE_FILES} ${dir}/${filename})
    ENDFOREACH(pathname)
ENDFOREACH(dir)

# collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)

FOREACH(subdir ${OPENMM_SOURCE_SUBDIRS})
    FILE(GLOB_RECURSE src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/${subdir}/src/*.c)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

SET(COMMON_KERNELS_CPP ${CMAKE_CURRENT_BINARY_DIR}/../common/src/CommonImageKernelSources.cpp)
SET(SOURCE_FILES ${SOURCE_FILES} ${COMMON_KERNELS_CPP})

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/../common/src)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/platforms/cuda/include)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/platforms/cuda/src)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_BINARY_DIR}/platforms/cuda/src)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/platforms/common/include)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_BINARY_DIR}/platforms/common/src)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR}/../common/src)

# Create the library

INCLUDE_DIRECTORIES(${CUDA_TOOLKIT_INCLUDE})
# GENERATED property. Is this source file generated as part of the build or CMake process.
SET_SOURCE_FILES_PROPERTIES(${COMMON_KERNELS_CPP} PROPERTIES GENERATED TRUE)
# there are class definations in .h files, so they should be parameters in the ADD_LIBRARY().
ADD_LIBRARY(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})
# 使顶级依赖于其他顶级<target>目标，以确保它们在构建之前<target>进行构建。
# 顶级目标是由 add_executable(), add_library()或add_custom_target()命令之一创建的(但不是CMake生成的目标，如install).
ADD_DEPENDENCIES(${SHARED_TARGET} ImageCommonKernels)

TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${CUDA_LIBRARIES})
TARGET_LINK_LIBRARIES(${SHARED_TARGET} OpenMM)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} OpenMMCUDA)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${IMAGE_LIBRARY_NAME})
# 设置目标的一些属性来改变它们构建的方式. COMPILE_FLAGS属性可以设置附加的编译器标志，它们会在构建目标内的源文件时被用到. -D是在预编译中定义一些变量
SET_TARGET_PROPERTIES(${SHARED_TARGET} PROPERTIES
    COMPILE_FLAGS "-DOPENMM_BUILDING_SHARED_LIBRARY ${EXTRA_COMPILE_FLAGS}"
    LINK_FLAGS "${EXTRA_COMPILE_FLAGS}")
IF (APPLE)
    SET_TARGET_PROPERTIES(${SHARED_TARGET} PROPERTIES LINK_FLAGS "-F/Library/Frameworks -framework CUDA ${EXTRA_COMPILE_FLAGS}")
ENDIF (APPLE)

INSTALL(TARGETS ${SHARED_TARGET} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/plugins)
# Ensure that links to the main CUDA library will be resolved.
IF (APPLE)
    SET(CUDA_LIBRARY libOpenMMCUDA.dylib)
    INSTALL(CODE "EXECUTE_PROCESS(COMMAND install_name_tool -change ${CUDA_LIBRARY} @loader_path/${CUDA_LIBRARY} ${CMAKE_INSTALL_PREFIX}/lib/plugins/lib${SHARED_TARGET}.dylib)")
ENDIF (APPLE)

add_subdirectory(tests)
