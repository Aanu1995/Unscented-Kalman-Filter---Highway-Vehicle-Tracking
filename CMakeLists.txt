cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

add_definitions(-std=c++17)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

project(playback)

find_package(PCL 1.15 REQUIRED)
find_package(Eigen3 3.4 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
link_directories(${PCL_LIBRARY_DIRS} ${EIGEN3_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS} ${EIGEN3_DEFINITIONS})
list(REMOVE_ITEM PCL_LIBRARIES "vtkproj4")


add_executable (ukf_highway src/main.cpp src/ukf.cpp src/tools.cpp src/render/render.cpp)
target_link_libraries (ukf_highway ${PCL_LIBRARIES} ${EIGEN3_LIBRARIES})




