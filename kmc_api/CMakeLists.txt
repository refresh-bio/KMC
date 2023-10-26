# -*- cmake -*-

project(KMC_API
        VERSION 3.2.2
        LANGUAGES C CXX)

# make list of all API source files
file(GLOB sources "[a-zA-Z]*.cpp")
      
add_library(kmc_api OBJECT
            ${sources})
