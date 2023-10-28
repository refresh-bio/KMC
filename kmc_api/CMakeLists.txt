# -*- cmake -*-

project(KMC_API
        VERSION 3.2.2
        LANGUAGES C CXX)

# make list of all API source files
file(GLOB sources "[a-zA-Z]*.cpp")
# export it to the parent scope, so python module could reuse it
set(PY_KMC_API_SOURCES ${sources} PARENT_SCOPE)
      
add_library(kmc_api OBJECT
            ${sources})
