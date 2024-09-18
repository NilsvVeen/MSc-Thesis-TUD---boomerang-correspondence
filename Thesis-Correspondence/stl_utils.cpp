
#include <iostream>
#include <filesystem>
#include <fstream>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <opencv2/opencv.hpp> // For saving the PNG

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>

#include <igl/readSTL.h>
#include <igl/writeSTL.h>
#include <igl/fit_plane.h>
#include <igl/slice_mask.h> // To filter the vertices



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

// Function to compute the surface area of a triangle
double computeTriangleArea(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) {
    return 0.5 * (v1 - v0).cross(v2 - v0).norm();
}


// Project vertices to 2D (XY plane) and fix Z to a specified value
Eigen::MatrixXd projectToXYWithFixedZ(const Eigen::MatrixXd& V, double fixedZ) {
    Eigen::MatrixXd V_2D(V.rows(), 3);
    V_2D.leftCols(2) = V.leftCols(2);  // Copy X and Y
    V_2D.col(2).setConstant(fixedZ);  // Set Z to fixedZ
    return V_2D;
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

    // Screenshot for the rotated mesh
    polyscope::screenshot("rotated_mesh.png");

    // Now we will slice the mesh along the Z-axis and calculate the area for each slice
    int numSlices = 10;
    double minZ = V.col(2).minCoeff();
    double maxZ = V.col(2).maxCoeff();
    double sliceThickness = (maxZ - minZ) / numSlices;

    int maxAreaSlice = -1;
    double maxSliceArea = 0.0;

    Eigen::MatrixXd V_maxAreaSlice;
    Eigen::MatrixXi F_maxAreaSlice;

    for (int i = 0; i < numSlices; ++i) {
        double sliceMinZ = minZ + i * sliceThickness;
        double sliceMaxZ = sliceMinZ + sliceThickness;

        std::vector<int> sliceFaces;  // Store the indices of faces in this slice

        // Identify faces that are within the slice range
        for (int j = 0; j < F.rows(); ++j) {
            Eigen::Vector3d v0 = V.row(F(j, 0));
            Eigen::Vector3d v1 = V.row(F(j, 1));
            Eigen::Vector3d v2 = V.row(F(j, 2));

            // Check if all vertices of this face are within the slice range
            if ((v0.z() >= sliceMinZ && v0.z() < sliceMaxZ) &&
                (v1.z() >= sliceMinZ && v1.z() < sliceMaxZ) &&
                (v2.z() >= sliceMinZ && v2.z() < sliceMaxZ)) {
                sliceFaces.push_back(j);  // Add the face index
            }
        }

        // Create the sliced faces matrix for this slice
        Eigen::MatrixXi F_slice(sliceFaces.size(), 3);
        for (int k = 0; k < sliceFaces.size(); ++k) {
            F_slice.row(k) = F.row(sliceFaces[k]);
        }

        // Compute the surface area for this slice
        double sliceArea = 0.0;
        for (int j = 0; j < F_slice.rows(); ++j) {
            Eigen::Vector3d v0 = V.row(F_slice(j, 0));
            Eigen::Vector3d v1 = V.row(F_slice(j, 1));
            Eigen::Vector3d v2 = V.row(F_slice(j, 2));
            sliceArea += computeTriangleArea(v0, v1, v2);
        }

        // Find the slice with the maximum surface area
        if (sliceArea > maxSliceArea) {
            maxSliceArea = sliceArea;
            maxAreaSlice = i;
            V_maxAreaSlice = V;
            F_maxAreaSlice = F_slice;
        }
    }

    if (maxAreaSlice >= 0) {
        // Compute the average Z coordinate for the max area slice
        Eigen::VectorXd zCoords(F_maxAreaSlice.rows());
        for (int i = 0; i < F_maxAreaSlice.rows(); ++i) {
            zCoords(i) = V.row(F_maxAreaSlice(i, 0)).z() +
                V.row(F_maxAreaSlice(i, 1)).z() +
                V.row(F_maxAreaSlice(i, 2)).z();
        }
        double avgZ = zCoords.mean() / 3.0;

        // Project the mesh onto the XY plane with fixed Z coordinate
        Eigen::MatrixXd V_2D = projectToXYWithFixedZ(V, avgZ);

        // Register the 2D projection with Polyscope
        polyscope::registerSurfaceMesh("2D Projection", V_2D, F);

        // Screenshot for the 2D projection
        polyscope::screenshot("2d_projection.png");
    }

    // Show the rotated mesh and 2D projection
    polyscope::show();
}