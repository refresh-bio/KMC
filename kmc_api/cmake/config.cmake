configure_file (cmake/<definitions.h>.in ${CMAKE_CURRENT_SOURCE_DIR}/<definitions.h>)

list (APPEND KMC_API_DEFINITIONS "my_fopen=fopen")
if (MSVC)
	list (APPEND KMC_API_DEFINITIONS
		"my_fseek=_fseeki64"
		"my_ftell=_ftelli64")
else (MSVC)
	list (APPEND KMC_API_DEFINITIONS
		"my_fseek=fseek"
		"my_ftell=ftell"
		"_TCHAR=char"
		"_tmain=main")
endif (MSVC)
# Version information
list (APPEND KMC_API_DEFINITIONS
	"KMC_VER=${PROJECT_VERSION_MAJOR}:${PROJECT_VERSION_MINOR}:${PROJECT_VERSION_PATCH}"
	"KMC_DATE=${PROJECT_VERSION_TWEAK}")

set_property (DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS "${KMC_API_DEFINITIONS}")
