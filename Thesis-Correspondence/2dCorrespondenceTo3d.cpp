#include <iostream>
#include <vector>
#include <tuple>
#include <Eigen/Core>
#include <filesystem>
#include "stl_utils.h"
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>


// Function to read meshes and point clouds
void readMeshesAndPointClouds(const std::string& meshesFolder,
    const std::string& pointCloudsFolder,
    Eigen::MatrixXd& Mesh1_V, Eigen::MatrixXi& Mesh1_F,
    Eigen::MatrixXd& Mesh2_V, Eigen::MatrixXi& Mesh2_F,
    std::vector<Eigen::MatrixXd>& V1_pointclouds,
    std::vector<Eigen::MatrixXd>& V2_pointclouds) {

    // Step 1: Read Meshes from the specified folder
    std::string mesh1File = meshesFolder + "/LeftMesh.obj";  // Example filenames
    std::string mesh2File = meshesFolder + "/RightMesh.obj";

    // Use utility functions to read meshes
    if (!readMeshFromFile(mesh1File, Mesh1_V, Mesh1_F)) {
        std::cerr << "Error reading mesh 1 from: " << mesh1File << std::endl;
        return;
    }

    if (!readMeshFromFile(mesh2File, Mesh2_V, Mesh2_F)) {
        std::cerr << "Error reading mesh 2 from: " << mesh2File << std::endl;
        return;
    }

    // Step 2: Read Pointclouds from the specified folder
    for (const auto& entry : std::filesystem::directory_iterator(pointCloudsFolder)) {
        std::string filePath = entry.path().string();

        if (filePath.find("V1__") != std::string::npos) {  // Read V1 point clouds
            Eigen::MatrixXd V1;
            if (readPointCloudFromFile(filePath, V1)) {
                V1_pointclouds.push_back(V1);
            }
            else {
                std::cerr << "Error reading V1 point cloud: " << filePath << std::endl;
            }
        }
        else if (filePath.find("V2_") != std::string::npos) {  // Read V2 point clouds
            Eigen::MatrixXd V2;
            if (readPointCloudFromFile(filePath, V2)) {
                V2_pointclouds.push_back(V2);
            }
            else {
                std::cerr << "Error reading V2 point cloud: " << filePath << std::endl;
            }
        }
    }
}

void showInPolyscope(const Eigen::MatrixXd& Mesh1_V, const Eigen::MatrixXi& Mesh1_F,
    const Eigen::MatrixXd& Mesh2_V, const Eigen::MatrixXi& Mesh2_F,
    const std::vector<Eigen::MatrixXd>& V1_pointclouds,
    const std::vector<Eigen::MatrixXd>& V2_pointclouds, bool enableMeshes = true, bool enablePointClouds = true) {

    // Initialize Polyscope
    polyscope::init();

    // Set color for the meshes
    glm::vec3 meshColor(1, 0, 0);

    // Register Mesh 1 (Surface Mesh) and set color
    auto mesh1_obj = polyscope::registerSurfaceMesh("Left Mesh", Mesh1_V, Mesh1_F);
    mesh1_obj->setSurfaceColor(meshColor);
    mesh1_obj->setEnabled(enableMeshes);

    // Register Mesh 2 (Surface Mesh) and set color
    auto mesh2_obj = polyscope::registerSurfaceMesh("Right Mesh", Mesh2_V, Mesh2_F);
    mesh2_obj->setSurfaceColor(meshColor);
    mesh2_obj->setEnabled(enableMeshes);

    // Register and configure point clouds, assigning different colors for each pair
    for (size_t i = 0; i < V1_pointclouds.size(); ++i) {
        // Generate a unique color for each i (using HSV to RGB or some other method)
        glm::vec3 pointCloudColor = glm::vec3(static_cast<float>(i) / V1_pointclouds.size(), 0.5f, 1.0f); // Example: color gradient

        // Register and configure V1 point cloud
        auto obj1 = polyscope::registerPointCloud("V1 Pointcloud " + std::to_string(i + 1), V1_pointclouds[i]);
        obj1->setPointColor(pointCloudColor);
        obj1->setPointRadius(0.005);
        obj1->setEnabled(enablePointClouds);

        // Register and configure V2 point cloud with the same color
        auto obj2 = polyscope::registerPointCloud("V2 Pointcloud " + std::to_string(i + 1), V2_pointclouds[i]);
        obj2->setPointColor(pointCloudColor);
        obj2->setPointRadius(0.005);
        obj2->setEnabled(enablePointClouds);
    }

    // Show Polyscope UI
    polyscope::show();
}

