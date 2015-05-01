set (KMC_COUNTER_EXTERNAL_INCLUDES
	${CMAKE_SOURCE_DIR}/external/asmlib.h
	${CMAKE_SOURCE_DIR}/external/bzlib.h
	${CMAKE_SOURCE_DIR}/external/bzlib_private.h
	${CMAKE_SOURCE_DIR}/external/zconf.h
	${CMAKE_SOURCE_DIR}/external/zlib.h)

if (WIN32 OR APPLE)
	set (KMC_COUNTER_EXTERNAL_LIBRARIES
		${CMAKE_SOURCE_DIR}/external/libz.1.2.5.dylib
		${CMAKE_SOURCE_DIR}/external/libbz2.1.0.5.dylib)
	if (${KMC_BUILD_WITH_ASM})
		list (APPEND KMC_COUNTER_EXTERNAL_LIBRARIES
			${CMAKE_SOURCE_DIR}/external/libamac64.a)
	endif()
else (WIN32 OR APPLE)
	set (KMC_COUNTER_EXTERNAL_LIBRARIES
		${CMAKE_SOURCE_DIR}/external/libz.a
		${CMAKE_SOURCE_DIR}/external/libbz2.a)
	if (${KMC_BUILD_WITH_ASM})
		list (APPEND KMC_COUNTER_EXTERNAL_LIBRARIES
			${CMAKE_SOURCE_DIR}/external/alibelf64.a)
	endif ()
endif (WIN32 OR APPLE)

if (${KMC_INSTALL_RUNTIME})
	add_executable (kmc_counter ${KMC_COUNTER_EXTERNAL_INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR}/kmc_counter.cpp)
	target_link_libraries (kmc_counter kmc_core ${KMC_COUNTER_EXTERNAL_LIBRARIES})
	install (TARGETS kmc_counter RUNTIME DESTINATION ${INSTALL_RUNTIME_DIR})
endif (${KMC_INSTALL_RUNTIME})
