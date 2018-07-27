cmake_minimum_required(VERSION 3.2.0)

project(retinal_dev)

find_package(BioDynaMo REQUIRED)
include(${BDM_USE_FILE})
include_directories("src")

file(GLOB_RECURSE HEADERS src/*.h)
file(GLOB_RECURSE SOURCES src/*.cc)

bdm_add_executable(retinal_dev
                   HEADERS ${HEADERS}
                   SOURCES ${SOURCES}
LIBRARIES ${BDM_REQUIRED_LIBRARIES})
