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


file(
	GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS
	${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_executable(
	S3D_test 
	${SRC_FILES}
	${CMAKE_CURRENT_SOURCE_DIR}/../sigma/lib/Specification/Specification.cpp
)

set(SGM_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../sigma/include)
set(EIGEN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../eigen)
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)

target_include_directories(
	S3D_test PRIVATE
	${SGM_DIR}
	${EIGEN_DIR}
	${INCLUDE_DIR}
)

target_compile_features(S3D_test PRIVATE cxx_std_17)
target_compile_definitions(S3D_test PRIVATE SGM_CXX_STANDARD=2017)