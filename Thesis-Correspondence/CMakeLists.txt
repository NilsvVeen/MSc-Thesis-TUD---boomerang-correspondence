# CMakeList.txt : CMake project for Thesis-Correspondence, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("Thesis-Correspondence")


# Add source to this project's executable.
add_executable (Thesis-Correspondence "Thesis-Correspondence.cpp" "Thesis-Correspondence.h"  "geometry.cpp" "geometry.h" "file_utils.cpp" "file_utils.h"   "stl_utils.cpp" "stl_utils.h" "border_vertex_extraction.cpp" "border_vertex_extraction.h" "parameterize.h" "parameterize.cpp" "2dCorrespondenceTo3d.h" "2dCorrespondenceTo3d.cpp" "parameterizeSurface.cpp" "parameterizeSurface.h" "random_files.cpp" "splitMesh.cpp" "splitMesh.h" "uvMapCorrespondence.cpp" "uvMapCorrespondence.h" "evaluateCorrespondence.h" "evaluateCorrespondence.cpp" "intermediateShape.h" "intermediateShape.cpp")

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Specify supported build types for multi-configuration systems (like Visual Studio)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Supported build configurations" FORCE)

# Optionally, add compiler flags or specific settings for each build type
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# Continue with the rest of the file...

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






# --- CGAL Setup ---
set(CGAL_DIR "${CMAKE_SOURCE_DIR}/libraries/CGAL-5.6.1")

# Manually set CMAKE_PREFIX_PATH for CGAL
set(CMAKE_PREFIX_PATH ${CGAL_DIR})

# Set the GMP directory
set(GMP_DIR "${CGAL_DIR}/auxiliary/gmp")
# Include GMP directories
set(GMP_INCLUDE_DIR "${GMP_DIR}/include")  # Adjust if necessary
set(GMP_LIBRARIES "${GMP_DIR}/lib/libgmp-10.lib")  # Adjust if necessary

# Set the MPFR directory
set(MPFR_DIR "${GMP_DIR}")  # Adjust if necessary, this should be where MPFR is located
set(MPFR_INCLUDE_DIR "${MPFR_DIR}/include")  # Adjust if necessary
set(MPFR_LIBRARIES "${MPFR_DIR}/lib/libmpfr-4.lib")  # Adjust if necessary

# Include directories for GMP and MPFR
include_directories(${GMP_INCLUDE_DIR} ${MPFR_INCLUDE_DIR})


# --- Boost Setup ---
set(BOOST_ROOT "${CMAKE_SOURCE_DIR}/libraries/boost_1_86_0_b1")
set(BOOST_LIBRARY_DIR "${BOOST_ROOT}/lib64-msvc-14.2")  # Adjust if necessary for your configuration
set(BOOST_INCLUDE_DIR "${BOOST_ROOT}")

# Find Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(Thesis-Correspondence ${Boost_LIBRARIES})


# Find CGAL
find_package(CGAL REQUIRED)

# Include CGAL directories
include_directories(${CGAL_INCLUDE_DIRS})

# Link CGAL, GMP, and MPFR libraries
target_link_libraries(Thesis-Correspondence ${CGAL_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES})




# # --- GLFW Setup ---
# include_directories("${CMAKE_SOURCE_DIR}/libraries/glfw-3.4")
# include_directories("${CMAKE_SOURCE_DIR}/libraries/glfw-3.4/include")
# target_link_libraries(Thesis-Correspondence glfw)


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
# target_link_libraries(Thesis-Correspondence igl::core)

# Enable libigl's CGAL functionality
set(LIBIGL_WITH_CGAL ON)
set(LIBIGL_WITH_CGAL_COPYLEFT ON)

# Include the CGAL copyleft directory
include_directories("${LIBIGL_DIR}/include/igl/copyleft/cgal")


# Link the necessary libigl components (e.g., core, OpenGL, GLFW, CGAL)
target_link_libraries(Thesis-Correspondence igl::core)


# --- OpenCV Setup ---

# Set the path to the directory containing OpenCVConfig.cmake
#set(OpenCV_DIR "${CMAKE_SOURCE_DIR}/libraries/opencv/build")

# Find OpenCV package
#find_package(OpenCV 4.10 REQUIRED)

# Add include directories for OpenCV
#include_directories(${OpenCV_INCLUDE_DIRS})


# Link OpenCV libraries to your target
#target_link_libraries(Thesis-Correspondence ${OpenCV_LIBS})

# Define a variable for the models directory
set(MODELS_DIRECTORY "${CMAKE_SOURCE_DIR}/models")
# Add the global models directory as a compile definition
add_definitions(-DMODELS_DIRECTORY="${MODELS_DIRECTORY}")