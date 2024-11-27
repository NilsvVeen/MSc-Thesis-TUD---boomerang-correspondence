#include <iostream>
#include <vector>
#include <tuple>
#include <Eigen/Core>
#include <filesystem>
#include "stl_utils.h"
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>

#include "stl_utils.h"
#include "file_utils.h"
#include <igl/point_mesh_squared_distance.h>


// Function to read meshes and point clouds
void readMeshesAndPointClouds(const std::string& meshesFolder,
    const std::string& pointCloudsFolder,
    Eigen::MatrixXd& Mesh1_V, Eigen::MatrixXi& Mesh1_F,
    Eigen::MatrixXd& Mesh2_V, Eigen::MatrixXi& Mesh2_F,
    std::vector<Eigen::MatrixXd>& V1_pointclouds,
    std::vector<Eigen::MatrixXd>& V2_pointclouds) {

    // Step 1: Read Meshes from the specified folder
    std::string mesh1File = meshesFolder + "/M1.obj";  // Example filenames
    std::string mesh2File = meshesFolder + "/M2.obj";

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


// Helper function to find exact (x, y) matches in Mesh1 and add the corresponding z value
void findExactCorrespondences(const Eigen::MatrixXd& mesh1_V, Eigen::MatrixXd& V1_pointcloud, std::vector<int>& indices_to_remove) {
    // Create a set to store already processed (x, y, z) points
    std::set<std::tuple<double, double, double>> processedPoints;

    // Temporary container for valid points
    Eigen::MatrixXd validPoints(V1_pointcloud.rows(), V1_pointcloud.cols());
    int validIndex = 0;

    
    for (int i = 0; i < V1_pointcloud.rows(); ++i) {
        bool matchFound = false;

        for (int j = 0; j < mesh1_V.rows(); ++j) {
            // Compare x, y coordinates
            if (V1_pointcloud(i, 0) == mesh1_V(j, 0) && V1_pointcloud(i, 1) == mesh1_V(j, 1)) {
                // Create the (x, y, z) tuple for the match
                std::tuple<double, double, double> point(
                    V1_pointcloud(i, 0),
                    V1_pointcloud(i, 1),
                    mesh1_V(j, 2)
                );

                // Check if the point is already added
                if (processedPoints.find(point) == processedPoints.end()) {
                    // Set z value from mesh1_V
                    V1_pointcloud(i, 2) = mesh1_V(j, 2);
                    processedPoints.insert(point);  // Mark this point as processed
                    matchFound = true;

                    // Add to valid points
                    validPoints.row(validIndex++) = V1_pointcloud.row(i);


                }
                break;  // Exit loop once a match is found
            }
        }
        if (!matchFound) {
            indices_to_remove.push_back(i);
            std::cout << "no match!?!?! at : " << std::to_string(V1_pointcloud(i)) << " at index " << i << std::endl;

        }

        //// Update progress
        //std::cout << "\rProcessing V1 point cloud: " << std::fixed << std::setprecision(2)
        //    << (static_cast<double>(i + 1) / V1_pointcloud.rows()) * 100 << "% completed" << std::flush;
    }

    // Resize validPoints to include only valid rows
    validPoints.conservativeResize(validIndex, V1_pointcloud.cols());
    V1_pointcloud = validPoints;  // Update the original point cloud

    std::cout << std::endl; // To move to the next line after the progress is complete
}


// Helper function to find the closest (x,y) match in Mesh2
void findClosestCorrespondences(const Eigen::MatrixXd& mesh2_V, Eigen::MatrixXd& V2_pointcloud) {
    for (int i = 0; i < V2_pointcloud.rows(); ++i) {
        double min_distance = std::numeric_limits<double>::max();
        int closest_idx = -1;

        for (int j = 0; j < mesh2_V.rows(); ++j) {
            // Calculate Euclidean distance in 2D (x, y)
            double distance = std::sqrt(std::pow(V2_pointcloud(i, 0) - mesh2_V(j, 0), 2) +
                std::pow(V2_pointcloud(i, 1) - mesh2_V(j, 1), 2));
            if (distance < min_distance) {
                min_distance = distance;
                closest_idx = j;
            }
        }

        // Use the closest point from Mesh2 to update the x, y, z values
        if (closest_idx != -1) {
            V2_pointcloud(i, 0) = mesh2_V(closest_idx, 0);
            V2_pointcloud(i, 1) = mesh2_V(closest_idx, 1);
            V2_pointcloud(i, 2) = mesh2_V(closest_idx, 2);
        }

        // Update progress
        std::cout << "\rProcessing V2 point cloud: " << std::fixed << std::setprecision(2)
            << (static_cast<double>(i + 1) / V2_pointcloud.rows()) * 100 << "% completed" << std::flush;
    }
    std::cout << std::endl; // To move to the next line after the progress is complete
}

// note works, but contains something wrong?
void projectAndSplitMesh(const Eigen::MatrixXd& mesh2_V, const Eigen::MatrixXi& mesh2_F,
    const Eigen::MatrixXd& V2_pointcloud, Eigen::MatrixXd& V2_pointcloud_new,
    Eigen::MatrixXd& updated_mesh_V, Eigen::MatrixXi& updated_mesh_F) {
    // Copy original mesh vertices and faces
    updated_mesh_V = mesh2_V;
    updated_mesh_F = mesh2_F;

    // Initialize V2_pointcloud_new as empty
    V2_pointcloud_new.resize(0, 3);

    for (int i = 0; i < V2_pointcloud.rows(); ++i) {
        // Find the closest points and their associated faces
        Eigen::VectorXd squared_distances;
        Eigen::MatrixXi closest_faces;
        Eigen::MatrixXd closest_vertices;
        igl::point_mesh_squared_distance(V2_pointcloud.row(i), mesh2_V, mesh2_F, squared_distances, closest_faces, closest_vertices);

        // Get the closest vertex on the mesh
        Eigen::RowVector3d projected_point = closest_vertices.row(0);  // closest vertex

        // Add the projected point as a new vertex
        int new_vertex_idx = updated_mesh_V.rows();
        updated_mesh_V.conservativeResize(updated_mesh_V.rows() + 1, Eigen::NoChange);
        updated_mesh_V.row(new_vertex_idx) = projected_point;

        // Add the projected point to V2_pointcloud_new
        V2_pointcloud_new.conservativeResize(V2_pointcloud_new.rows() + 1, Eigen::NoChange);
        V2_pointcloud_new.row(V2_pointcloud_new.rows() - 1) = projected_point;

        // Find the face that the point was closest to
        int face_idx = closest_faces(0, 0);  // Assuming the point is closest to a single face

        // Get the vertices of the closest face
        Eigen::RowVector3i face = updated_mesh_F.row(face_idx);

        // Split the face by adding the new vertex
        Eigen::RowVector3i new_face1(face(0), face(1), new_vertex_idx);
        Eigen::RowVector3i new_face2(face(1), face(2), new_vertex_idx);
        Eigen::RowVector3i new_face3(face(2), face(0), new_vertex_idx);

        // Remove the old face and add the new faces
        updated_mesh_F.row(face_idx) = updated_mesh_F.row(updated_mesh_F.rows() - 1);  // Replace with last face
        updated_mesh_F.conservativeResize(updated_mesh_F.rows() - 1, Eigen::NoChange);  // Remove last face
        updated_mesh_F.conservativeResize(updated_mesh_F.rows() + 3, Eigen::NoChange);  // Add 3 new faces
        updated_mesh_F.row(updated_mesh_F.rows() - 3) = new_face1;
        updated_mesh_F.row(updated_mesh_F.rows() - 2) = new_face2;
        updated_mesh_F.row(updated_mesh_F.rows() - 1) = new_face3;

        // Update progress
        std::cout << "\rProcessing V2 point cloud: " << std::fixed << std::setprecision(2)
            << (static_cast<double>(i + 1) / V2_pointcloud.rows()) * 100 << "% completed" << std::flush;
    }

    std::cout << std::endl; // To move to the next line after the progress is complete
}

void projectAndReplaceVertices(const Eigen::MatrixXd& mesh2_V, const Eigen::MatrixXi& mesh2_F,
    const Eigen::MatrixXd& V2_pointcloud, Eigen::MatrixXd& V2_pointcloud_new,
    Eigen::MatrixXd& updated_mesh_V, Eigen::MatrixXi& updated_mesh_F) {
    // Copy original mesh vertices and faces
    updated_mesh_V = mesh2_V;
    updated_mesh_F = mesh2_F;

    // Initialize V2_pointcloud_new as empty
    V2_pointcloud_new.resize(0, 3);

    for (int i = 0; i < V2_pointcloud.rows(); ++i) {
        // Find the closest points and their associated faces
        Eigen::VectorXd squared_distances;
        Eigen::MatrixXi closest_faces;
        Eigen::MatrixXd closest_vertices;
        igl::point_mesh_squared_distance(V2_pointcloud.row(i), mesh2_V, mesh2_F, squared_distances, closest_faces, closest_vertices);

        // Get the projected point on the mesh
        Eigen::RowVector3d projected_point = closest_vertices.row(0);  // Closest point on mesh

        // Find the closest vertex in the mesh to the projected point
        int closest_vertex_idx = -1;
        double min_distance = std::numeric_limits<double>::max();

        for (int j = 0; j < updated_mesh_V.rows(); ++j) {
            double distance = (updated_mesh_V.row(j) - projected_point).squaredNorm();
            if (distance < min_distance) {
                min_distance = distance;
                closest_vertex_idx = j;
            }
        }

        if (closest_vertex_idx == -1) {
            std::cerr << "Error: No closest vertex found for projected point at index " << i << std::endl;
            continue;
        }

        // Replace the closest vertex with the projected point
        updated_mesh_V.row(closest_vertex_idx) = projected_point;

        // Add the projected point to V2_pointcloud_new
        V2_pointcloud_new.conservativeResize(V2_pointcloud_new.rows() + 1, Eigen::NoChange);
        V2_pointcloud_new.row(V2_pointcloud_new.rows() - 1) = projected_point;

        // Update progress
        std::cout << "\rProcessing V2 point cloud: " << std::fixed << std::setprecision(2)
            << (static_cast<double>(i + 1) / V2_pointcloud.rows()) * 100 << "% completed" << std::flush;
    }

    std::cout << std::endl; // To move to the next line after the progress is complete
}



// Function to write outputs to a folder
void writeOutputsToFolder(const std::string& outputFolder,
    const Eigen::MatrixXd& Mesh1_V, const Eigen::MatrixXi& Mesh1_F,
    const Eigen::MatrixXd& Mesh2_V, const Eigen::MatrixXi& Mesh2_F,
    const std::vector<Eigen::MatrixXd>& V1_pointclouds,
    const std::vector<Eigen::MatrixXd>& V2_pointclouds) {

    // Create the output folder if it doesn't exist
    createDirectory(outputFolder);

    // Optionally clear the directory first
    clearDirectory(outputFolder);

    // Define output filenames for the meshes
    std::string mesh1Filename = outputFolder + "/Mesh1.obj";
    std::string mesh2Filename = outputFolder + "/Mesh2.obj";

    // Save Mesh1 and Mesh2 using the existing save function
    saveMeshToFile(mesh1Filename, Mesh1_V, Mesh1_F);
    saveMeshToFile(mesh2Filename, Mesh2_V, Mesh2_F);

    // Write point clouds to files
    for (size_t i = 0; i < V1_pointclouds.size(); ++i) {
        std::string V1_filename = outputFolder + "/V1_pointcloud_" + std::to_string(i) + ".txt";
        savePointCloudToFile(V1_filename, V1_pointclouds[i]);  // Use your existing function

        std::string V2_filename = outputFolder + "/V2_pointcloud_" + std::to_string(i) + ".txt";
        savePointCloudToFile(V2_filename, V2_pointclouds[i]);  // Use your existing function
    }

    std::cout << "Outputs saved to: " << outputFolder << std::endl;
}
