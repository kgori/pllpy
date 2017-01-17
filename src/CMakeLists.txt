cmake_minimum_required(VERSION 3.2)
project(src)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")

set(SOURCE_FILES
    pllml.cpp
    pllml.h
    test.cpp)

set(MY_LIB_LINK_LIBRARIES -lpll-sse3)
add_executable(chk ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(chk ${MY_LIB_LINK_LIBRARIES})
