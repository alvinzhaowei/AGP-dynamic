cmake_minimum_required(VERSION 3.25)
project(DAGP)

set(CMAKE_CXX_STANDARD 17)

add_executable(DAGP main.cpp graph/Vertex.cpp graph/Graph.cpp  MyLib/ParaReader.cpp MyLib/MyTimer.cpp
        DS/DS.cpp
        DS/DS.h)