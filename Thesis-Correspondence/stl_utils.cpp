
#include <iostream>
#include <filesystem>
#include <igl/readSTL.h>
#include <fstream>
#include <Eigen/Core>
#include <polyscope/polyscope.h>
#include <opencv2/opencv.hpp> // For saving the PNG
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <igl/writeSTL.h>
#include <igl/fit_plane.h>


void viewSTLObject(const std::string& filename) {
    Eigen::MatrixXd V; // Vertices
    Eigen::MatrixXi F; // Faces
    Eigen::MatrixXd N; // Normals (if required)

    // Open the STL file using ifstream
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Failed to open STL file: " << filename << std::endl;
        return;
    }

    // Read the STL file using libigl
    if (!igl::readSTL(file, V, F, N)) {
        std::cerr << "Failed to load STL file: " << filename << std::endl;
        return;
    }

    // Close the file
    file.close();

    // Initialize the Polyscope viewer
    polyscope::init();

    // Register the mesh to Polyscope
    polyscope::registerSurfaceMesh("STL Mesh", V, F);

    // Show the viewer
    polyscope::show();
}

void printFilesInDirectory(const std::string& directoryPath) {
    try {
        std::filesystem::path dir(directoryPath);
        if (!std::filesystem::exists(dir)) {
            std::cerr << directoryPath << " does not exist." << std::endl;
            return;
        }
        if (!std::filesystem::is_directory(dir)) {
            std::cerr << directoryPath << " is not a valid directory." << std::endl;
            return;
        }

        std::cout << "Files in directory " << directoryPath << ":" << std::endl;
        for (const auto& entry : std::filesystem::directory_iterator(dir)) {
            std::cout << entry.path().filename().string() << std::endl;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error accessing directory: " << e.what() << std::endl;
    }
}


void fitPlaneAndAlignMesh(const std::string& filename) {
    Eigen::MatrixXd V; // Vertices
    Eigen::MatrixXi F; // Faces
    Eigen::MatrixXd N; // Normals

    // Open the STL file using ifstream
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Failed to open STL file: " << filename << std::endl;
        return;
    }

    // Read the STL file using libigl
    if (!igl::readSTL(file, V, F, N)) {
        std::cerr << "Failed to load STL file: " << filename << std::endl;
        return;
    }

    // Fit a plane to the mesh
    Eigen::RowVector3d planeNormal;
    Eigen::RowVector3d planePoint;
    igl::fit_plane(V, planeNormal, planePoint);

    // Output the plane normal
    std::cout << "Plane Normal: " << planeNormal << std::endl;

    // Calculate the rotation matrix to align planeNormal with the targetDirection (Z-axis)
    Eigen::Vector3d normal(planeNormal);
    Eigen::Vector3d targetDirection(0, 0, 1);
    Eigen::Matrix3d rotationMatrix = Eigen::Quaterniond::FromTwoVectors(normal, targetDirection).toRotationMatrix();

    // Rotate the mesh
    V = (rotationMatrix * V.transpose()).transpose();

    // Initialize Polyscope
    polyscope::init();

    // Register the rotated mesh with Polyscope
    polyscope::registerSurfaceMesh("Rotated Mesh", V, F);

    // Show the rotated mesh
    polyscope::show();
}