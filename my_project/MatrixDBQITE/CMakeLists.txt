cmake_minimum_required(VERSION 3.20)

project(MatrixDBQITE)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)

add_executable(dbqite
    main.cpp
    functions.cpp
    data.cpp

)

target_link_libraries(dbqite Eigen3::Eigen)
