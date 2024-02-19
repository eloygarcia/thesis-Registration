cmake_minimum_required(VERSION 2.8)
project(MRImammoRegistration)


OPTION(USE_CUDA "Build GPU Ray-cast (requires CUDA)" OFF)

IF( USE_CUDA )

# ------------------------------------------------------------
if(COMMAND cmake_policy)
   cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

IF (MSVC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_SCL_SECURE_NO_WARNINGS")
ENDIF (MSVC)

IF (NOT DEFINED CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Default build type, optimisation enabled.")
ENDIF (NOT DEFINED CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/source)

# -------------------------------------------------------------

find_package(ITK)
include(${ITK_USE_FILE})

find_package(ItkVtkGlue)
include(${ItkVtkGlue_USE_FILE})

# -------------------------------------------------------------
find_package(VTK)
include(${VTK_USE_FILE})

get_directory_property(dir_defs DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS)
set(vtk_flags)
foreach(it ${dir_defs})
    if(it MATCHES "vtk*")
    list(APPEND vtk_flags ${it})
    endif()
endforeach()

foreach(d ${vtk_flags})
    remove_definitions(-D${d})
endforeach()

# -------------------------------------------------------------
find_package(CUDA REQUIRED)
include_directories(
						${CUDA_INCLUDE_DIRS}
						${CUDA_SDK_COMMON_INCLUDE_DIR}
						)

cuda_include_directories(	
							${CUDA_INCLUDE_DIRS}	
							${CUDA_SDK_COMMON_INCLUDE_DIR}
							)

find_path(CUDA_SDK_COMMON_INCLUDE_DIR
      helper_cuda.h
      PATHS ${CUDA_SDK_SEARCH_PATH}
      PATH_SUFFIXES "common/inc"
      DOC "Location of helper_cuda.h"
     NO_DEFAULT_PATH
      )


CUDA_INCLUDE_DIRECTORIES(
							${CMAKE_SOURCE_DIR}/src/
							${CMAKE_SOURCE_DIR}/source
							)

cuda_compile(HELPER_O helper.cu)
cuda_compile(DRIVER_O driver.cu OPTIONS -G)

set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ;-D_GPU_" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GPU_" )
#  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} --gpu-architecture sm_20)
# set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} --gpu-architecture sm_20 )
set(CUDA_SEPARABLE_COMPILATION OFF)  # En off evita un error de compilación del cuda !

set(CUDA_VERBOSE_BUILD ON CACHE BOOL "nvcc verbose" )
include(FindCUDA)

# SDK ROOT => c:\Program Data\Nvidia Corporation\Cuda Sample \v.6.0.
# Toolkit Root => C:\Program Files\Nvidia GPU Computing toolkit\CUDA\v.6.0

# -------------------------------------------------------------
IF (WIN32 AND NOT CYGWIN)
  SET(DEF_INSTALL_CMAKE_DIR CMake)
ELSE ()
  SET(DEF_INSTALL_CMAKE_DIR lib/CMake/CudaProj)
ENDIF (WIN32 AND NOT CYGWIN)

# -------------------------------------------------------------
cuda_add_library(	xRayMammograms
					mammogram.cpp	)

# cuda_add_library(	MagneticResonanceImages
#					MRimage.cpp
#					changeSegmentation.cpp
#					mesh.cpp
#					cleaningMesh.cpp
#					tetgen.cxx
#					predicates.cxx	)

add_library(	Optimization
				Transformations.cpp
				IntensityBased.cpp
#				HillClimbingOpt.cpp	
#				StochasticHillClimbing.cpp
				SimulatedAnnealing.cpp
				Metrics.cpp
				MechanicalProperties.cpp)

add_library(	NiftySimModels
				NiftySimEjecutable.cpp
				xml.cpp
				xmlModelWriter.cpp	)

cuda_add_library(	RayCastingProjections
# 					classicProjection.cpp
					RegularGrid.cpp
#					NewProjection.cpp
					CudaProjection.cu
					CudaLibrary_kernels.cu	)

cuda_add_executable(	MRImammoRegistration
						# auxiliarDefinitions.cpp
						Source.cpp	)


target_link_libraries(	MRImammoRegistration
						${ITK_LIBRARIES} 
						${VTK_LIBRARIES} 
						${Glue}
						xRayMammograms
						# MagneticResonanceImages
						NiftySimModels
						RayCastingProjections
						Optimization
						)

# -------------------------------------------------------------
ELSE()  # USE_CUDA

find_package(ITK)
include(${ITK_USE_FILE})

find_package(ItkVtkGlue)
include(${ItkVtkGlue_USE_FILE})

find_package(VTK)
include(${VTK_USE_FILE})

add_library(	xRayMammograms
					mammogram.cpp	)

# add_library(	MagneticResonanceImages
#					MRimage.cpp
#					changeSegmentation.cpp
#					mesh.cpp
#					cleaningMesh.cpp
#					tetgen.cxx
#					predicates.cxx	)

add_library(	Optimization
				Transformations.cpp
				IntensityBased.cpp
#				HillClimbingOpt.cpp	
#				StochasticHillClimbing.cpp
				SimulatedAnnealing.cpp 
				Metrics.cpp
				MechanicalProperties.cpp)

add_library(	NiftySimModels
				NiftySimEjecutable.cpp
				xml.cpp
				xmlModelWriter.cpp	)

add_library(	RayCastingProjections
					RegularGrid.cpp
					# classicProjection.cpp
					# NewProjection.cpp	)
					)

add_executable(	MRImammoRegistration
						# auxiliarDefinitions.cpp
						Source.cpp	)


target_link_libraries(	MRImammoRegistration
						${ITK_LIBRARIES} 
						${VTK_LIBRARIES} 
						${Glue}
						xRayMammograms
						# MagneticResonanceImages
						NiftySimModels
						RayCastingProjections
						Optimization
						)

ENDIF( USE_CUDA )