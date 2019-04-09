cmake_minimum_required(VERSION 3.0)
project(WSA)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
set(CMAKE_CXX_FLAGS "-Wall -msse2 -g3 -O3 -funroll-all-loops -mavx2 -I ./ -std=c++11 -pthread")
set(CMAKE_C_FLAGS "-Wall -msse2 -g3 -O3 -funroll-all-loops -mavx2 -I ./ -pthread")

#for google unit test
add_subdirectory(./googletest/googletest)
include_directories(googletest/googletest/include googletest/googletest)


add_executable(wsa main.cpp
        src/seed.cpp src/seed.h
        src/fasta.cpp src/fasta.h
        src/test/seed_test.cpp
        src/zdp.cpp src/zdp.h src/Score.cpp src/Score.h src/gffToCategory.cpp src/gffToCategory.h src/SuffixTree.cpp src/SuffixTree.h src/test/suffixTree_test.cpp)

target_link_libraries(wsa gtest gtest_main)
