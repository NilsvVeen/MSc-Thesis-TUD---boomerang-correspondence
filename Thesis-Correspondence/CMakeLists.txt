# CMakeList.txt : CMake project for Thesis-Correspondence, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("Thesis-Correspondence")

# Add source to this project's executable.
add_executable (Thesis-Correspondence "Thesis-Correspondence.cpp" "Thesis-Correspondence.h"  "geometry.cpp" "geometry.h" "file_utils.cpp" "file_utils.h" "sfml_utils.cpp" "sfml_utils.h" "stl_utils.cpp" "stl_utils.h")

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


# --- Polyscope Setup ---
add_subdirectory("${CMAKE_SOURCE_DIR}/libraries/polyscope")
include_directories("${CMAKE_SOURCE_DIR}/libraries/polyscope/include")
target_link_libraries(Thesis-Correspondence polyscope)

# --- libigl Setup ---
set(LIBIGL_DIR "${CMAKE_SOURCE_DIR}/libraries/libigl")
add_subdirectory(${LIBIGL_DIR} "libigl")

# Ensure Eigen is included (comes with libigl)
include_directories("${LIBIGL_DIR}/include")

# Enable libigl OpenGL and GLFW components
set(LIBIGL_USE_STATIC_LIBRARY ON)        # Enable static linking (optional)
set(LIBIGL_WITH_OPENGL ON)               # Enable OpenGL
set(LIBIGL_WITH_OPENGL_GLFW ON)          # Enable GLFW for windowing
set(LIBIGL_WITH_OPENGL_GLFW_IMGUI ON)    # Enable ImGui for GUI interaction (optional)

# Link the necessary libigl components (e.g., core, OpenGL, and GLFW)
target_link_libraries(Thesis-Correspondence igl::core)



# Define a variable for the models directory
set(MODELS_DIRECTORY "${CMAKE_SOURCE_DIR}/models")

# Add the global models directory as a compile definition
add_definitions(-DMODELS_DIRECTORY="${MODELS_DIRECTORY}")