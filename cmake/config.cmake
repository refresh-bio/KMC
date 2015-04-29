configure_file (${PROJECT_SOURCE_DIR}/cmake/definitions.h.in definitions.h)

### Debug
list (APPEND KMC_DEBUG_DEFINITIONS
		"DEBUG_MODE"
		"memcpy"
		"memset")
set_property (DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS_DEBUG "${KMC_DEBUG_DEFINITIONS}")

### Develop
list (APPEND KMC_DEVELOP_DEFINITIONS "DEVELOP_MODE")
set_property (DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS_DEVELOP "${KMC_DEBUG_DEFINITIONS}")

### 
if (${KMC_DISABLE_ASM_LIB})
	list (APPEND KMC_DEFINITIONS "DISABLE_ASMLIB")
endif (${KMC_DISABLE_ASM_LIB})

### Required
list (APPEND KMC_DEFINITIONS "my_fopen=fopen")
if (MSVC)
	list (APPEND KMC_DEFINITIONS
		"my_fseek=_fseeki64"
		"my_ftell=_ftelli64")
else (MSVC)
	list (APPEND KMC_DEFINITIONS
		"my_fseek=fseek"
		"my_ftell=ftell"
		"_TCHAR=char"
		"_tmain=main")
endif (MSVC)
# Version information
list (APPEND KMC_DEFINITIONS
	"KMC_VER=${PROJECT_VERSION_MAJOR}${PROJECT_VERSION_MINOR}${PROJECT_VERSION_PATCH}"
	"KMC_DATE=${PROJECT_VERSION_TWEAK}")
# Section needs heading
list (APPEND KMC_DEFINITIONS
	"_CRT_SECURE_NO_WARNINGS"
	"USE_META_PROG"
	"KMER_X=3"
	"MAX_K=256"
	"MIN_K=10"
	"MIN_MEM=1")
# Range of number of bins
list (APPEND KMC_DEFINITIONS
	"MIN_N_BINS=64"
	"MAX_N_BINS=2000")
# Range of number of FASTQ/FASTA reading threads
list (APPEND KMC_DEFINITIONS
	"MIN_SF=1"
	"MAX_SF=32")
# Range of number of signature length
list (APPEND KMC_DEFINITIONS
	"MIN_SL=5"
	"MAX_SL=8")
# Range of number of splitting threads
list (APPEND KMC_DEFINITIONS
	"MIN_SP=1"
	"MAX_SP=64")
# Range of number of sorting threads
list (APPEND KMC_DEFINITIONS
	"MIN_SO=1"
	"MAX_SO=64")
# Range of number of sorter threads pre sorter in strict memory mode
list (APPEND KMC_DEFINITIONS
	"MIN_SMSO=1"
	"MAX_SMSO=16")
# Range of number of uncompactor threads in strict memory mode
list (APPEND KMC_DEFINITIONS
	"MIN_SMUN=1"
	"MAX_SMUN=16")
# Range of number of merger threads in strict memory mode
list (APPEND KMC_DEFINITIONS
	"MIN_SMME=1"
	"MAX_SMME=16")
# Range of number of threads per single sorting thread
list (APPEND KMC_DEFINITIONS
	"MIN_SR=1"
	"MAX_SR=16")

set_property (DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS "${KMC_DEFINITIONS}")
