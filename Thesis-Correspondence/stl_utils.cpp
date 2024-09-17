
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



// Function to compute the weighted average normal direction, rotate the object, and show it
void computeAndSaveNormalDirection(const std::string& filename) {
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

    // Compute the weighted average normal direction
    Eigen::Vector3d averageNormal(0, 0, 0);
    double totalArea = 0.0;

    for (int i = 0; i < F.rows(); ++i) {
        Eigen::Vector3d v0 = V.row(F(i, 0));
        Eigen::Vector3d v1 = V.row(F(i, 1));
        Eigen::Vector3d v2 = V.row(F(i, 2));

        // Compute the normal for this face
        Eigen::Vector3d normal = N.row(i).normalized();

        // Compute the area of the triangle
        Eigen::Vector3d e1 = v1 - v0;
        Eigen::Vector3d e2 = v2 - v0;
        double area = (e1.cross(e2)).norm() / 2.0;

        // Weight the normal by the area of the triangle
        totalArea += area;
        averageNormal += normal * area;
    }

    if (totalArea > 0) {
        averageNormal /= totalArea;
        averageNormal.normalize();
    }
    else {
        std::cerr << "Total area is zero. Cannot compute average normal direction." << std::endl;
        return;
    }

    std::cout << "Weighted Average Normal Direction: " << averageNormal.transpose() << std::endl;

    // Calculate the rotation matrix to align averageNormal with the targetDirection (Z-axis)
    Eigen::Vector3d targetDirection(0, 0, 1);
    Eigen::Matrix3d rotationMatrix = Eigen::Quaterniond::FromTwoVectors(averageNormal, targetDirection).toRotationMatrix();

    // Rotate the mesh
    V = (rotationMatrix * V.transpose()).transpose();

    // Initialize Polyscope
    polyscope::init();

    // Register the rotated mesh with Polyscope
    polyscope::registerSurfaceMesh("Rotated Mesh", V, F);

    // Show the rotated mesh
    polyscope::show();
}