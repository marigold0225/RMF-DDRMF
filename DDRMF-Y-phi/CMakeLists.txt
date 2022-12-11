cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(recipe-09 LANGUAGES Fortran)
include(GNUInstallDirs)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
# defines targets and sources
#add_subdirectory(src)
#add_subdirectory(lib)
# contains an "external" library we will link to
# enable testing and define tests
#function(LinkFortranLibraries Target)
#	target_link_libraries(Target ${ARGN})
#	foreach(Lib IN LISTS ARGN)
#		target_include_directories(run PUBLIC $<TARGET_PROPERTY:${Lib},Fortran_MODULE_DIRECTORY>)
#	endforeach()
#endfunction()a
#这当前项目里lib目录下的所有文件(包含子目录下的文件)添加到SRCS这个变量里
file(GLOB_RECURSE SRCS RELATIVE ${PROJECT_SOURCE_DIR} lib/*)
# 将SRCS编译为mylib这个动态链接库
add_library(mylib SHARED ${SRCS})

file(GLOB_RECURSE SRCS RELATIVE ${PROJECT_SOURCE_DIR} src/*)
add_executable(run ${SRCS})
# 让myprogram依赖mylib
add_dependencies(run mylib)
# 链接mylib
target_link_libraries(run mylib)
