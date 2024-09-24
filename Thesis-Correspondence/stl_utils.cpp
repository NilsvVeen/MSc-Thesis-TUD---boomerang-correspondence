#include <iostream>
#include <filesystem>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Dense>

//#include <opencv2/opencv.hpp> // For saving the PNG

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>

#include <igl/readSTL.h>
#include <igl/writeSTL.h>
#include <igl/fit_plane.h>
#include <igl/slice_mask.h> // To filter the vertices
#include <igl/copyleft/cgal/convex_hull.h>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <vector>
#include <polyscope/pick.h>

#include <algorithm> // For std::find



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


// Save mesh to a file in OBJ format
void saveMeshToFile(const std::string& filename, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    // Write vertices
    for (int i = 0; i < V.rows(); ++i) {
        file << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
    }

    // Write faces
    for (int i = 0; i < F.rows(); ++i) {
        file << "f " << F(i, 0) + 1 << " " << F(i, 1) + 1 << " " << F(i, 2) + 1 << std::endl;
    }

    file.close();
    std::cout << "Mesh saved to " << filename << std::endl;
}

// Save point cloud to a file in PLY format
void savePointCloudToFile(const std::string& filename, const Eigen::MatrixXd& V) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    // Write PLY header
    file << "ply\nformat ascii 1.0\nelement vertex " << V.rows() << "\n";
    file << "property float x\nproperty float y\nproperty float z\nend_header\n";

    // Write vertices
    for (int i = 0; i < V.rows(); ++i) {
        file << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
    }

    file.close();
    std::cout << "Point cloud saved to " << filename << std::endl;
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



bool readMeshFromFile(const std::string& filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file for reading: " << filename << std::endl;
        return false;
    }

    std::vector<Eigen::Vector3d> vertices;
    std::vector<Eigen::Vector3i> faces;
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string type;
        iss >> type;

        if (type == "v") {
            Eigen::Vector3d vertex;
            iss >> vertex[0] >> vertex[1] >> vertex[2];
            vertices.push_back(vertex);
        }
        else if (type == "f") {
            Eigen::Vector3i face;
            int idx;
            for (int i = 0; i < 3; ++i) {
                iss >> idx;
                face[i] = idx - 1;  // OBJ format uses 1-based indexing
            }
            faces.push_back(face);
        }
    }

    V = Eigen::MatrixXd(vertices.size(), 3);
    for (size_t i = 0; i < vertices.size(); ++i) {
        V.row(i) = vertices[i];
    }

    F = Eigen::MatrixXi(faces.size(), 3);
    for (size_t i = 0; i < faces.size(); ++i) {
        F.row(i) = faces[i];
    }

    file.close();
    return true;
}





bool readPointCloudFromFile(const std::string& filename, Eigen::MatrixXd& V) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file for reading: " << filename << std::endl;
        return false;
    }

    std::string line;
    bool headerEnded = false;
    std::vector<Eigen::Vector3d> vertices;

    while (std::getline(file, line)) {
        if (!headerEnded) {
            if (line == "end_header") {
                headerEnded = true;
            }
            continue;
        }

        std::istringstream iss(line);
        Eigen::Vector3d vertex;
        iss >> vertex[0] >> vertex[1] >> vertex[2];
        vertices.push_back(vertex);
    }

    V = Eigen::MatrixXd(vertices.size(), 3);
    for (size_t i = 0; i < vertices.size(); ++i) {
        V.row(i) = vertices[i];
    }

    file.close();
    return true;
}


void showPointCloudInParts(const Eigen::MatrixXd& V) {

    int splits = 10;

    int totalVertices = V.rows();
    int partSize = totalVertices / splits; // 20% of the point cloud

    polyscope::init();

    // Loop through and register each part of the point cloud
    for (int i = 0; i < splits; ++i) {
        Eigen::MatrixXd partV = V.middleRows(i * partSize, partSize);
        std::string cloudName = "Point Cloud Part " + std::to_string(i + 1);
        polyscope::registerPointCloud(cloudName, partV);
    }

    polyscope::show();
}


void showSelection(const Eigen::MatrixXd& V) {
    // Initialize Polyscope
    polyscope::init();

    // Create a vector to keep track of selected vertex indices
    std::vector<int> selectedVertices;

    // Register the point cloud with Polyscope
    auto* pointCloud = polyscope::registerPointCloud("Point Cloud", V);

    // Enable picking and handle selection updates in the callback
    polyscope::state::userCallback = [&]() {
        // Check if there is a current selection
        if (polyscope::pick::haveSelection()) {
            // Get the selection pair (object, localIndex)
            auto pickResult = polyscope::pick::getSelection();

            // Check if the selected object is the point cloud
            if (pickResult.first == pointCloud) {
                int selectedIndex = pickResult.second; // Index of the selected vertex

                // Check if the vertex is already in the list
                auto it = std::find(selectedVertices.begin(), selectedVertices.end(), selectedIndex);
                if (it != selectedVertices.end()) {
                    // If found, remove it (deselect)
                    selectedVertices.erase(it);

                    // Optional: Print the deselected vertex index
                    std::cout << "Deselected vertex: " << selectedIndex << std::endl;
                }
                else {
                    // If not found, add it (select)
                    selectedVertices.push_back(selectedIndex);

                    // Optional: Print the selected vertex index
                    std::cout << "Selected vertex: " << selectedIndex << std::endl;
                }
            }
        }

        // Optional: Print all selected vertices
        std::cout << "Currently selected vertices: ";
        for (int v : selectedVertices) {
            std::cout << v << " ";
        }
        std::cout << std::endl;
    };

    // Show the Polyscope UI
    polyscope::show();
}