cmake_minimum_required(VERSION 3.16)

project(LTI_project VERSION 0.1.0)

#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
#if(COMPILER_SUPPORTS_MARCH_NATIVE)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
#endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")
elseif (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  /std:c++17")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++17 -O3")
endif (CMAKE_COMPILER_IS_GNUCXX)


set(PROJECT_TOP_DIR ${PROJECT_SOURCE_DIR})
SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/bin)

FILE(GLOB SRC "${PROJECT_SOURCE_DIR}/src/*.cpp")
include_directories(${PROJECT_SOURCE_DIR}/Include)
include_directories(${PROJECT_SOURCE_DIR}/lib/eigen)

add_library(LTI INTERFACE)
add_executable(develop main.cpp)

option(ENABLE_UNIT_TESTS "Enable unit tests" ON)
option(ENABLE_LCOV "Enable code coverage tests" OFF)
message(STATUS "---------------------------")
message(STATUS "Current : ${PROJECT_NAME}")
message(STATUS "Enable testing: ${ENABLE_UNIT_TESTS}")
message(STATUS "Enable coverage: ${ENABLE_LCOV}")

if(${ENABLE_UNIT_TESTS})
    add_definitions(-DENABLE_UNIT_TESTS)
    add_subdirectory(test)
endif()


# copy all lib to Twincat SDK
if(DEFINED ENV{TWINCATSDK})
    file(COPY Include DESTINATION $ENV{TWINCATSDK})
    MESSAGE(STATUS "Copy to Twincat SDK: $ENV{TWINCATSDK}")
endif()