#	SPDX-FileCopyrightText: (c) 2025 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
#	SPDX-License-Identifier: MIT License
#=========#=========#=========#=========#=========#=========#=========#=========#=========#=========


include(FetchContent)
FetchContent_Declare(
	How2use
	GIT_REPOSITORY https://github.com/SigmaPark/How2use.git
	GIT_TAG v2.1.1
)

set(BUILD_HOW2USE_TEST_PROJ OFF CACHE BOOL "Disable How2use's test project build" FORCE)
FetchContent_MakeAvailable(How2use)

if(MSVC)
	add_definitions(-D_WINDOWS_SOLUTION_)

	add_compile_options(/bigobj)
	add_definitions(-DUNICODE)
	add_definitions(-D_UNICODE)

endif()

if(UNIX)
	add_definitions(-D_POSIX_SOLUTION_)

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

message("-- CMAKE_C_COMPILER=${CMAKE_C_COMPILER}")
message("-- CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

file(
	GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS
	${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/*.hpp
)

add_executable(S3D_test ${SRC_FILES})

target_link_libraries(S3D_test PRIVATE S3D_lib)
target_link_libraries(S3D_test PRIVATE How2use_lib)