#	SPDX-FileCopyrightText: (c) 2021 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
#	SPDX-License-Identifier: MIT License
#=========#=========#=========#=========#=========#=========#=========#=========#=========#=========


include(ExternalProject)
set(EIGEN_BUILD_DIR ${CMAKE_BINARY_DIR}/eigen)
set(EIGEN_LIB_DIR ${EIGEN_BUILD_DIR}/src/eigen)
ExternalProject_Add(
	eigen
	URL https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
	PREFIX ${EIGEN_BUILD_DIR}
	CONFIGURE_COMMAND ""
	BUILD_COMMAND ""
	INSTALL_COMMAND ""
	TEST_COMMAND ""
)

include(FetchContent)
FetchContent_Declare(
	sigma
	GIT_REPOSITORY https://github.com/SigmaPark/homo-genericus.git
	GIT_TAG v2.0.0
)

set(BUILD_SIGMA_TEST_PROJ OFF CACHE BOOL "Disable sigma's test project build" FORCE)
FetchContent_MakeAvailable(sigma)

set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)

file(
	GLOB_RECURSE INCLUDE_FILES CONFIGURE_DEPENDS
	${INCLUDE_DIR}/*.hpp
)

add_library(S3D_lib INTERFACE ${INCLUDE_FILES})
add_dependencies(S3D_lib eigen)

target_link_libraries(S3D_lib INTERFACE Sigma_lib)

target_include_directories(
	S3D_lib INTERFACE 
	${INCLUDE_DIR}
	${EIGEN_LIB_DIR}
)