#	SPDX-FileCopyrightText: (c) 2021 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
#	SPDX-License-Identifier: MIT License
#=========#=========#=========#=========#=========#=========#=========#=========#=========#=========


set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)
set(SGM_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../sigma/include)
set(EIGEN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../eigen)

file(
	GLOB_RECURSE INCLUDE_FILES CONFIGURE_DEPENDS
	${INCLUDE_DIR}/*.hpp
)

add_library(S3D_lib INTERFACE ${INCLUDE_FILES})

target_include_directories(
    S3D_lib INTERFACE 
    ${INCLUDE_DIR}
	${SGM_INCLUDE_DIR}
	${EIGEN_DIR}    
)

target_compile_features(S3D_lib INTERFACE cxx_std_17)
target_compile_definitions(S3D_lib INTERFACE SGM_CXX_STANDARD=2017)