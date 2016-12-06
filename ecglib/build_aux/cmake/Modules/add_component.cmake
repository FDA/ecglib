function(add_component componentname dependentcomponents)

project(${componentname})

set(INCS "")
set(LIBS "")
foreach(comp ${dependentcomponents})
	if(EXISTS ${${comp}_SOURCE_DIR})
		set(INCS ${INCS} ${${comp}_SOURCE_DIR})
	else()
		message(WARNING "${comp}_SOURCE_DIR not defined!")
	endif()

	if(TARGET ${lib${comp}})
		set(LIBS ${LIBS} ${lib${comp}})
	else()
		message(WARNING "lib${comp} is not a target!")
	endif()
endforeach(comp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${INCS})

file(GLOB_RECURSE CODE_FILES *.cpp)

IF( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" AND NOT CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
	# Base
	install(DIRECTORY . DESTINATION include/ FILES_MATCHING PATTERN "*.hpp")
	install(DIRECTORY . DESTINATION include/ FILES_MATCHING PATTERN "*.h")
ENDIF()

if(UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()

if(WITH_STATIC)
	add_library(${PROJECT_NAME}-static STATIC ${CODE_FILES})

	IF( CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
		set_target_properties(${PROJECT_NAME}-static PROPERTIES OUTPUT_NAME ${PROJECT_NAME}-dbg)
	ELSE( CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
		set_target_properties(${PROJECT_NAME}-static PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
	ENDIF()

	target_link_libraries(${PROJECT_NAME}-static ${ECGLIB_LIBRARIES} ${LIBS})
	
	install(TARGETS ${PROJECT_NAME}-static ARCHIVE DESTINATION lib)
endif()

if(WITH_SHARED)
	add_library(${PROJECT_NAME} SHARED ${CODE_FILES})

	IF( CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
		set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME}-dbg)
	ELSE( CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
		set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
	ENDIF()

	target_link_libraries(${PROJECT_NAME} ${ECGLIB_LIBRARIES} ${LIBS})

	if(WIN32)
		install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION lib)
	else()
		install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION lib)
	endif()
endif()

set(lib${PROJECT_NAME} ${PROJECT_NAME} CACHE INTERNAL "lib${PROJECT_NAME}")
if(NOT WITH_SHARED)
	set(lib${PROJECT_NAME} ${PROJECT_NAME}-static CACHE INTERNAL "lib${PROJECT_NAME}")
endif()

endfunction(add_component)
