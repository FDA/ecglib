# # ecglib - C++ library for processing and analysis of electrocardiograms (ECGs)
#
# ## LICENSE
#
# Authors - Lars Johannesen, Jose Vicente, Meisam Hosseini and Dustin C McAfee
# ecglib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License (GNU LGPL3), or (at your option) any later version.
# ecglib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see http://www.gnu.org/licenses/.
#
# ## DISCLAIMER
# ecglib software and documentation were developed by the authors in their capacities as  Oak Ridge Institute for Science and Education (ORISE) research fellows at the U.S. Food and Drug Administration (FDA). .
# FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic.  Further, FDA makes no representations that the use of the Software will not infringe any patent or proprietary rights of third parties.   The use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions.

cmake_minimum_required(VERSION 2.8)

########################################################################################
include_directories(src/core/)
########################################################################################
# Detect compiler flags for C++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -DARMA_DONT_USE_CXX11")
else()
        message(ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

CHECK_CXX_COMPILER_FLAG("-ftrack-macro-expansion=0" COMPILER_NEEDS_DISABLEMACROEXPANSION)
if(COMPILER_NEEDS_DISABLEMACROEXPANSION)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftrack-macro-expansion=0")
else()
	message("Compiler is ok with macros")
endif()

########################################################################################
# Determine version
set(VER_MAJOR 1)
set(VER_MINOR 0)
set(VER_PATCH 0)

########################################################################################
# Determine build type
option(WITH_STATIC "With static" OFF)
if(NOT WITH_STATIC)
	if(UNIX)
		set(WITH_STATIC 0)
	else()
		set(WITH_STATIC 1)
	endif()
endif()

option(WITH_SHARED "With shared" OFF)
if(NOT WITH_SHARED)
	if(UNIX)
		set(WITH_SHARED 1)
	else()
		set(WITH_SHARED 0)
	endif()
endif()

if(NOT WITH_SHARED AND NOT WITH_STATIC)
	message(FATAL_ERROR "Cannot compile without static and shared!")
endif()

option(WITH_TEST "With testing" OFF)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "Setting build type to 'Release' as none was specified.")
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
	# Set the possible values of build type for cmake-gui
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

########################################################################################
# Detect necessary libraries
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build_aux/cmake/Modules/")
set(ECGLIB_LIBRARIES "" CACHE INTERNAL "ECGLIB LIBRARIES")

########################################################################################
# Configure options
set(TYPES CORE;DELINEATORS)
set(COMPONENTS-CORE core)
set(COMPONENTS-DELINEATORS delineators-twave)

# Build options
set(ALL_COMPONENTS "")
set(SELECTED_COMPONENTS "")
foreach(typ ${TYPES})
	set(WITH_COMPONENTS_${typ} ${COMPONENTS-${typ}} CACHE STRING "${typ} components")
	set(ALL_COMPONENTS ${ALL_COMPONENTS} ${COMPONENTS-${typ}})
	set(SELECTED_COMPONENTS ${SELECTED_COMPONENTS} ${WITH_COMPONENTS_${typ}})
endforeach(typ)

# Order them appropriately
set(BUILD_COMPONENTS "")

# For all possible components
foreach(acomp ${ALL_COMPONENTS})
	# Is this one selected
	foreach(scomp ${SELECTED_COMPONENTS})
		# If yes, add it as the next component
		if(${acomp} STREQUAL ${scomp})
			set(BUILD_COMPONENTS ${BUILD_COMPONENTS} ${scomp})
		endif()
	endforeach(scomp)
endforeach(acomp)

##################################################
# Boost
if( "$ENV{BOOST_ROOT}" STREQUAL "")
	message(STATUS "environment variable BOOST_ROOT not defined searching in default paths")
else()
	message(STATUS "BOOST_ROOT set to: $ENV{BOOST_ROOT}, so skipping system path")
	set(Boost_NO_SYSTEM_PATHS TRUE)
endif()

include(FindBoost)
find_package(Boost COMPONENTS date_time filesystem system regex program_options iostreams serialization)

if(Boost_FOUND)
	message(STATUS "Boost_MAJOR_VERSION = ${Boost_MAJOR_VERSION}")
	message(STATUS "Boost_MINOR_VERSION = ${Boost_MINOR_VERSION}")
	message(STATUS "Boost_INCLUDE_DIR = ${Boost_INCLUDE_DIR}")
	message(STATUS "Boost_LIB_DIR = ${Boost_LIBRARY_DIRS}")

	# FIXME Check version has fs3
	include_directories(include ${Boost_INCLUDE_DIR})
	set(ECGLIB_LIBRARIES "${ECGLIB_LIBRARIES};${Boost_LIBRARIES}" CACHE INTERNAL "ECGLIB LIBRARIES")
else()
	message(FATAL_ERROR "Boost not found")
endif()

##################################################
# Armadillo (lapack/blas)
include(FindArmadillo)
if(ARMADILLO_FOUND)
	message(STATUS "Armadillo_MAJOR_VERSION = ${ARMADILLO_VERSION_MAJOR}")
	message(STATUS "Armadillo_MINOR_VERSION = ${ARMADILLO_VERSION_MINOR}")

	include_directories(include ${ARMADILLO_INCLUDE_DIRS})
	set(ECGLIB_LIBRARIES "${ECGLIB_LIBRARIES};${ARMADILLO_LIBRARY}" CACHE INTERNAL "ECGLIB LIBRARIES")
else()
	message(FATAL_ERROR "Armadillo not found")
endif()

if( "$ENV{BLAS_LAPACK}" STREQUAL "")
	include(FindBLAS)
	include(FindLAPACK)
	if(BLAS_FOUND AND LAPACK_FOUND)
		message(STATUS "Blas libraries = ${BLAS_LIBRARIES}")
		message(STATUS "Lapack libraries = ${LAPACK_LIBRARIES}")

		set(ECGLIB_LIBRARIES "${ECGLIB_LIBRARIES};${LAPACK_LIBRARIES}" CACHE INTERNAL "ECGLIB LIBRARIES")
	else()
		message(FATAL_ERROR "Lapack/Blas not found")
	endif()
else()
	message(STATUS "Skipping autodetect for lapack/atlas and use the following libraries: ${$ENV{LAPACK_ROOT}}")
	set(ECGLIB_LIBRARIES "${ECGLIB_LIBRARIES};$ENV{BLAS_LAPACK}" CACHE INTERNAL "ECGLIB LIBRARIES")
endif()

if(WIN32)
	set(ECGLIB_LIBRARIES "${ECGLIB_LIBRARIES};-lgfortran;-lblas;-llapack;-lquadmath;-lws2_32")
endif()

########################################################################################
# Setup packaging / output
IF( CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
	SET(PACKAGE_NAME "ecglib-dbg")
	SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -gdwarf-3 ${CMAKE_CXX_DEBUG_FLAGS} -DTIXML_USE_TICPP=YES")
ELSE( CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
	SET(PACKAGE_NAME "ecglib")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTIXML_USE_TICPP=YES")
ENDIF()

set(CPACK_PACKAGE_DESCRIPTION "ECG Processing Library")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ECG Processing Library written in C++")
set(CPACK_PACKAGE_NAME ${PACKAGE_NAME})
set(CPACK_DEBIAN_PACKAGE_DEPENDS "")
set(CPACK_PACKAGE_CONTACT "Jose Vicente")
set(CPACK_PACKAGE_VENDOR "U.S. Food and Drug Administration")
set(CPACK_PACKAGE_VERSION_MAJOR ${VER_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${VER_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${VER_PATCH})


########################################################################################
# Generate config.hpp
if(WIN32)
	set(CPACK_GENERATOR "TGZ")
	set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
else()
	set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}-${CMAKE_SYSTEM_PROCESSOR}")
	set(CPACK_GENERATOR "DEB")
endif()

include(CPack)

message(STATUS "**************************SUMMARY*************************")
message(STATUS "Version: ${ECGLIB_VERSION}")
message(STATUS "Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Testing: ${WITH_TEST}")
message(STATUS "Static: ${WITH_STATIC}")
message(STATUS "Shared: ${WITH_SHARED}")
message(STATUS "With:")

# Build them
set(DOC_DIRS_LIST "")
set(DOC_OVERVIEW_DIRS_LIST "")
foreach(bcomp ${BUILD_COMPONENTS})
	string(REGEX REPLACE "-" "/" cdir ${bcomp})
	set(cdir "src/${cdir}")
	add_subdirectory(${cdir})
	message(STATUS "\t${bcomp}")

	set(DOC_DIRS_LIST ${DOC_DIRS_LIST} "${CMAKE_SOURCE_DIR}/${cdir}")

	set(ddir "${CMAKE_SOURCE_DIR}/${cdir}/doc")
	if(EXISTS ${ddir})
		set(DOC_DIRS_LIST ${DOC_DIRS_LIST} ${ddir})

		set(ddir "${CMAKE_SOURCE_DIR}/${cdir}/doc/overview.dox")
		if(EXISTS ${ddir})
			set(ddir "@ref ${bcomp}")
			set(DOC_OVERVIEW_DIRS_LIST ${DOC_OVERVIEW_DIRS_LIST} ${ddir})
		endif()
	endif()
endforeach(bcomp)

string(REPLACE ";" " " DOC_DIRS "${DOC_DIRS_LIST}")
string(REPLACE ";" "<br>" DOC_OVERVIEW_DIRS "${DOC_OVERVIEW_DIRS_LIST}")

message(STATUS "**********************************************************")

########################################################################################
# Docs
find_package(Doxygen)
if(DOXYGEN_FOUND)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/ecglib.dox.in ${CMAKE_CURRENT_SOURCE_DIR}/doc/ecglib.dox @ONLY)
	add_custom_target(doc
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM
	)
endif(DOXYGEN_FOUND)
