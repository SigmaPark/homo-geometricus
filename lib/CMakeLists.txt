#	SPDX-FileCopyrightText: (c) 2021 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
#	SPDX-License-Identifier: MIT License
#=========#=========#=========#=========#=========#=========#=========#=========#=========#=========


set(EXTERNAL_LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../external_library)

include(ExternalProject)

set(EIGEN_BUILD_DIR   ${CMAKE_BINARY_DIR}/eigen)
set(EIGEN_INSTALL_DIR ${EXTERNAL_LIBRARY_DIR}/Eigen)
set(EIGEN_INCLUDE_DIR ${EIGEN_INSTALL_DIR})

ExternalProject_Add(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
    PREFIX ${EIGEN_BUILD_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND
      ${CMAKE_COMMAND} -E copy_directory ${EIGEN_BUILD_DIR}/src/eigen/Eigen ${EIGEN_INCLUDE_DIR}
    TEST_COMMAND ""
)
include_directories(SYSTEM ${EIGEN_INCLUDE_DIR})

set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)
set(SGM_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../sigma/include)

file(
	GLOB_RECURSE INCLUDE_FILES CONFIGURE_DEPENDS
	${INCLUDE_DIR}/*.hpp
)

add_library(S3D_lib INTERFACE ${INCLUDE_FILES})

target_include_directories(
    S3D_lib INTERFACE 
    ${INCLUDE_DIR}
	${SGM_INCLUDE_DIR}	
	${EXTERNAL_LIBRARY_DIR}
)

target_compile_features(S3D_lib INTERFACE cxx_std_17)
target_compile_definitions(S3D_lib INTERFACE SGM_CXX_STANDARD=2017)