# CMakeList.txt : CMake project for Thesis-Correspondence, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("Thesis-Correspondence")

# Add source to this project's executable.
add_executable (Thesis-Correspondence "Thesis-Correspondence.cpp" "Thesis-Correspondence.h" "curve_parametrization.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Thesis-Correspondence PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
