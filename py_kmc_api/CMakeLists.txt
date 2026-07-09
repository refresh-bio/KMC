# -*- cmake -*-

project(PY_KMC_API
        LANGUAGES C CXX)

find_package(Python COMPONENTS Interpreter Development REQUIRED)

# FIXME: switch to pybind11 cmake
Python_add_library(py_kmc_api MODULE
  ${PY_KMC_API_SOURCES}
  py_kmc_api.cpp)
target_include_directories(py_kmc_api PRIVATE "libs/pybind11/include")
# FIXME: a bit hacky
target_include_directories(py_kmc_api PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../kmc_api)

install(TARGETS py_kmc_api
        DESTINATION bin
        COMPONENT python_module)
install(FILES py_kmc_dump.py
        DESTINATION bin
        COMPONENT python_module)
