cmake_minimum_required(VERSION 3.0)
project(SeebeckL CXX)

# Fail if cmake is called in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

set(CMAKE_CXX_STANDARD 11)
add_compile_options("-O3")

# add 3rd-party libraries (before setting compiler flags etc)
include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/cubature")
add_subdirectory(3rdparty)

find_package(GSL REQUIRED)
if (GSL_FOUND)
    include_directories(GSL_INCLUDE_DIRS)
    link_libraries(GSL_LIBRARIES)
else()
    message(FATAL_ERROR "Not Found GSL Library, Please install GSL or manually \
        specify the installation directory for GSL, use \"-DGSL_ROOT_DIR=[...]\"")
endif()

# 添加头文件目录
include_directories("./include")