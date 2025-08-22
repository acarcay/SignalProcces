cmake_minimum_required(VERSION 3.15)
project(SignalProcces)

set(CMAKE_CXX_STANDARD 17)

find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})

add_executable(SignalProcces main.cpp filters.cpp)

target_link_libraries(SignalProcces PRIVATE ${Python3_LIBRARIES})
