add_subdirectory (api)
add_subdirectory (core)

if (${KMC_BUILD_RUNTIME})
	add_subdirectory (count)
	add_subdirectory (dump)
	add_subdirectory (sample)
endif (${KMC_BUILD_RUNTIME})
