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

# Specify SFML directory
set(SFML_DIR "${CMAKE_SOURCE_DIR}/libraries/SFML-2.6.1/lib/cmake/SFML")

# Find SFML package
find_package(SFML 2.6 REQUIRED COMPONENTS graphics window system)

# Include directories for SFML
include_directories("${CMAKE_SOURCE_DIR}/libraries/SFML-2.6.1/include")

# Link SFML libraries
if (SFML_FOUND)
    target_link_libraries(Thesis-Correspondence 
        sfml-graphics 
        sfml-window 
        sfml-system
    )
else()
    message(FATAL_ERROR "SFML not found. Please check your SFML installation and path.")
endif()