cmake_minimum_required(VERSION 3.17)
project(SoftProj_Final C)

set(CMAKE_C_STANDARD 90)
find_package(PythonLibs REQUIRED)


find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

add_executable(SoftProj_Final spkmeans.c )

add_library(spkmeans.h spkmeansmodule.c spkmeans.c)

target_link_libraries(spkmeans.h spkmeansmodule.c ${PYTHON_LIBRARIES})