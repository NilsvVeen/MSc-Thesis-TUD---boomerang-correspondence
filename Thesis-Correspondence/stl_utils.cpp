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

// Initialize Polyscope and register a mesh
void initializePolyscopeAndRegisterMesh(const std::string& name, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    polyscope::init();
    polyscope::registerSurfaceMesh(name, V, F);
}

// Take a screenshot of the current Polyscope view
void takeScreenshot(const std::string& screenshotName) {
    polyscope::screenshot(screenshotName);
}

// Show all registered meshes in Polyscope
void showMeshes() {
    polyscope::show();
}

// Function to fit a plane to the mesh and align it
void fitPlaneAndAlignMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& rotatedV, Eigen::MatrixXi& rotatedF) {
    Eigen::RowVector3d planeNormal;
    Eigen::RowVector3d planePoint;
    igl::fit_plane(V, planeNormal, planePoint);

    Eigen::Vector3d normal(planeNormal);
    Eigen::Vector3d targetDirection(0, 0, 1);
    Eigen::Matrix3d rotationMatrix = Eigen::Quaterniond::FromTwoVectors(normal, targetDirection).toRotationMatrix();

    rotatedV = (rotationMatrix * V.transpose()).transpose();
    rotatedF = F;
}

// Function to find the slice with maximum area
void findMaxAreaSlice(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& V_maxAreaSlice, Eigen::MatrixXi& F_maxAreaSlice, int& maxAreaSlice) {
    int numSlices = 10;
    double minZ = V.col(2).minCoeff();
    double maxZ = V.col(2).maxCoeff();
    double sliceThickness = (maxZ - minZ) / numSlices;

    double maxSliceArea = 0.0;
    maxAreaSlice = -1;

    for (int i = 0; i < numSlices; ++i) {
        double sliceMinZ = minZ + i * sliceThickness;
        double sliceMaxZ = sliceMinZ + sliceThickness;

        std::vector<int> sliceFaces;

        for (int j = 0; j < F.rows(); ++j) {
            Eigen::Vector3d v0 = V.row(F(j, 0));
            Eigen::Vector3d v1 = V.row(F(j, 1));
            Eigen::Vector3d v2 = V.row(F(j, 2));

            if ((v0.z() >= sliceMinZ && v0.z() < sliceMaxZ) &&
                (v1.z() >= sliceMinZ && v1.z() < sliceMaxZ) &&
                (v2.z() >= sliceMinZ && v2.z() < sliceMaxZ)) {
                sliceFaces.push_back(j);
            }
        }

        Eigen::MatrixXi F_slice(sliceFaces.size(), 3);
        for (int k = 0; k < sliceFaces.size(); ++k) {
            F_slice.row(k) = F.row(sliceFaces[k]);
        }

        double sliceArea = 0.0;
        for (int j = 0; j < F_slice.rows(); ++j) {
            Eigen::Vector3d v0 = V.row(F_slice(j, 0));
            Eigen::Vector3d v1 = V.row(F_slice(j, 1));
            Eigen::Vector3d v2 = V.row(F_slice(j, 2));
            sliceArea += computeTriangleArea(v0, v1, v2);
        }

        if (sliceArea > maxSliceArea) {
            maxSliceArea = sliceArea;
            maxAreaSlice = i;
            V_maxAreaSlice = V;
            F_maxAreaSlice = F_slice;
        }
    }
}

// Function to calculate average Z coordinate for a slice
double calculateAverageZ(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    Eigen::VectorXd zCoords(F.rows());
    for (int i = 0; i < F.rows(); ++i) {
        zCoords(i) = V.row(F(i, 0)).z() +
            V.row(F(i, 1)).z() +
            V.row(F(i, 2)).z();
    }
    return zCoords.mean() / 3.0;
}

// Function to handle STL file loading and processing
void processSTLFile(const std::string& filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& N) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Failed to open STL file: " << filename << std::endl;
        return;
    }

    if (!igl::readSTL(file, V, F, N)) {
        std::cerr << "Failed to load STL file: " << filename << std::endl;
        return;
    }

    file.close();
}

// Main function to fit plane, align mesh, and show results
void fitPlaneAndAlignMesh(const std::string& filename) {
    Eigen::MatrixXd V; // Vertices
    Eigen::MatrixXi F; // Faces
    Eigen::MatrixXd N; // Normals

    // Process STL file
    processSTLFile(filename, V, F, N);

    // Fit plane and align mesh
    Eigen::MatrixXd rotatedV;
    Eigen::MatrixXi rotatedF;
    fitPlaneAndAlignMesh(V, F, rotatedV, rotatedF);

    // Initialize Polyscope and register the rotated mesh
    initializePolyscopeAndRegisterMesh("Rotated Mesh", rotatedV, rotatedF);
    takeScreenshot("rotated_mesh.png");

    // Find the slice with maximum area
    Eigen::MatrixXd V_maxAreaSlice;
    Eigen::MatrixXi F_maxAreaSlice;
    int maxAreaSlice;
    findMaxAreaSlice(rotatedV, rotatedF, V_maxAreaSlice, F_maxAreaSlice, maxAreaSlice);

    if (maxAreaSlice >= 0) {
        // Calculate average Z coordinate for the max area slice
        double avgZ = calculateAverageZ(rotatedV, F_maxAreaSlice);

        // Project the mesh onto the XY plane with fixed Z coordinate
        Eigen::MatrixXd V_2D = projectToXYWithFixedZ(rotatedV, avgZ);

        // Register the 2D projection and max area slice with Polyscope
        initializePolyscopeAndRegisterMesh("2D Projection", V_2D, rotatedF);
        initializePolyscopeAndRegisterMesh("Max Area Slice", V_maxAreaSlice, F_maxAreaSlice);

        // Take screenshots
        takeScreenshot("2d_projection.png");
        takeScreenshot("max_area_slice.png");
    }

    // Show the registered meshes
    showMeshes();
}


