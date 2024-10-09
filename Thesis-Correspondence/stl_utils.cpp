#include <iostream>
#include <filesystem>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <tuple>

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
#include <polyscope/curve_network.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "parameterize.h"




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


// unused
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
    pointCloud->setPointRadius(0.0005); // Set the radius to 0.005

    // Create an array for default vertex colors (white)
    std::vector<std::array<double, 3>> vertexColors(V.rows(), { {1.0, 1.0, 1.0} }); // Default color: white

    // Add vertex color quantity to the point cloud
    pointCloud->addColorQuantity("vertex color", vertexColors);

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
                    polyscope::pick::resetSelection();
                }
                else {
                    // If not found, add it (select)
                    selectedVertices.push_back(selectedIndex);
                    polyscope::pick::resetSelection();

                }

                // Update vertex colors based on selection
                std::vector<std::array<double, 3>> updatedColors = vertexColors; // Copy current colors
                for (int idx : selectedVertices) {
                    updatedColors[idx] = { {0.0, 1.0, 0.0} }; // Set selected vertex color to green
                }

                // Apply the updated colors to the point cloud
                pointCloud->addColorQuantity("vertex color", updatedColors);
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


// unused
void showSideBySideSelection(const Eigen::MatrixXd& V1, const Eigen::MatrixXd& V2) {
    // Initialize Polyscope
    polyscope::init();

    // Register the first point cloud
    auto* pointCloud1 = polyscope::registerPointCloud("Point Cloud 1", V1);
    pointCloud1->setPointRadius(0.005); // Set point radius

    // Compute the bounding box for V1
    Eigen::VectorXd minV1 = V1.colwise().minCoeff();
    Eigen::VectorXd maxV1 = V1.colwise().maxCoeff();
    double widthV1 = maxV1(0) - minV1(0);  // Width along the X-axis

    // Compute the bounding box for V2
    Eigen::VectorXd minV2 = V2.colwise().minCoeff();
    Eigen::VectorXd maxV2 = V2.colwise().maxCoeff();
    double widthV2 = maxV2(0) - minV2(0);  // Width along the X-axis

    // Calculate the average width of the two point clouds
    double avgWidth = (widthV1 + widthV2) / 2.0;

    // Add additional space based on the average width divided by 10
    double extraSpace = avgWidth / 10.0;

    // Set the offset to ensure no overlap, adding extra space
    double offset = widthV1 + extraSpace;

    // Offset the second point cloud by translating it along the X-axis
    Eigen::MatrixXd V2_offset = V2;
    V2_offset.col(0).array() += offset;

    // Register the second point cloud with the offset applied
    auto* pointCloud2 = polyscope::registerPointCloud("Point Cloud 2", V2_offset);
    pointCloud2->setPointRadius(0.005); // Set point radius for the second point cloud

    // Show the Polyscope UI
    polyscope::show();
}





void updateVertexColors(polyscope::PointCloud* pointCloud, const std::vector<int>& selectedVertices, const std::vector<std::array<double, 3>>& defaultColors) {
    // Create a copy of the current colors
    std::vector<std::array<double, 3>> updatedColors = defaultColors;

    // Set the selected vertices to green
    for (int idx : selectedVertices) {
        updatedColors[idx] = { {0.0, 1.0, 0.0} }; // Green for selected vertices
    }

    // Update the color quantity in Polyscope
    pointCloud->addColorQuantity("vertex color", updatedColors);
}

void handleSelection(polyscope::PointCloud* pointCloud, std::vector<int>& selectedVertices, const std::vector<std::array<double, 3>>& defaultColors) {
    int selectedIndex = polyscope::pick::getSelection().second;

    // Check if the vertex is already selected
    auto it = std::find(selectedVertices.begin(), selectedVertices.end(), selectedIndex);
    if (it != selectedVertices.end()) {
        // Deselect if it's already selected
        selectedVertices.erase(it);
    }
    else {
        // Otherwise, select the vertex
        selectedVertices.push_back(selectedIndex);
    }

    // Reset the selection and update the colors
    polyscope::pick::resetSelection();
    updateVertexColors(pointCloud, selectedVertices, defaultColors);
}

polyscope::PointCloud* registerPointCloudWithColors(const std::string& name, const Eigen::MatrixXd& V, double pointRadius, std::vector<std::array<double, 3>>& vertexColors) {
    // Register the point cloud
    auto* pointCloud = polyscope::registerPointCloud(name, V);
    pointCloud->setPointRadius(pointRadius);

    // Set the initial colors (all white)
    vertexColors = std::vector<std::array<double, 3>>(V.rows(), { {1.0, 1.0, 1.0} });
    pointCloud->addColorQuantity("vertex color", vertexColors);

    return pointCloud;
}

Eigen::MatrixXd offsetPointCloud(const Eigen::MatrixXd& V, double offset) {
    Eigen::MatrixXd V_offset = V;
    V_offset.col(0).array() += offset;
    return V_offset;
}


void updateLineConnections(const Eigen::MatrixXd& V1, const Eigen::MatrixXd& V2_offset, const std::vector<int>& selectedVertices1, const std::vector<int>& selectedVertices2, const double& radius_default) {
    // Determine the number of matching pairs based on the smaller size
    int numLines = std::min(selectedVertices1.size(), selectedVertices2.size());

    if (numLines > 0) {
        std::vector<std::array<int, 2>> edges;
        Eigen::MatrixXd points(numLines * 2, 3); // Each line has 2 points

        // Create line segments for matching points in both point clouds
        for (int i = 0; i < numLines; ++i) {
            int idx1 = selectedVertices1[i];
            int idx2 = selectedVertices2[i];

            // Line goes from V1[idx1] to V2_offset[idx2]
            points.row(2 * i) = V1.row(idx1);
            points.row(2 * i + 1) = V2_offset.row(idx2);

            // Add line segment between consecutive points
            edges.push_back({ 2 * i, 2 * i + 1 });
        }

        // Register or update the curve network (lines) between the selected vertices
        auto* curveNetwork = polyscope::registerCurveNetwork("Corresponding Lines", points, edges);

        // Set the radius of the curve network
        curveNetwork->setRadius(radius_default);
        
    }
}


Eigen::MatrixXd calculateAndAdjustOffsets(const Eigen::MatrixXd& V1, const Eigen::MatrixXd& V2) {
    // Calculate necessary offsets
    double maxX_V1 = V1.col(0).maxCoeff(); // Maximum X of V1
    double minX_V2 = V2.col(0).minCoeff(); // Minimum X of V2

    // Calculate the offset to ensure V2 starts after V1
    double offsetX = maxX_V1 - minX_V2;

    // Offset V2 to be next to V1, keeping the z-coordinates the same
    Eigen::MatrixXd V2_offset = V2;

    // Adjust V2's X coordinates to be next to V1
    V2_offset.col(0) = V2.col(0).array() + offsetX; // Shift in x-direction

    // Keep z-coordinates the same for alignment
    if (V2_offset.rows() <= V1.rows()) {
        V2_offset.col(2) = V1.col(2).head(V2_offset.rows()); // Ensure z-coordinates are the same
    }

    // Align the Y-coordinates
    double averageY1 = V1.col(1).mean();
    V2_offset.col(1) = V2.col(1).array() + (averageY1 - V2.col(1).mean()); // Center V2 around the average Y of V1

    return V2_offset;
}

void handleUserSelection(auto pointCloud1, auto pointCloud2,
    const Eigen::MatrixXd& V1, const Eigen::MatrixXd& V2_offset,
    std::vector<int>& selectedVertices1, std::vector<int>& selectedVertices2,
    const std::vector<std::array<double, 3>>& vertexColors1,
    const std::vector<std::array<double, 3>>& vertexColors2,
    double radius_default) {
    if (polyscope::pick::haveSelection()) {
        auto pickResult = polyscope::pick::getSelection();

        if (pickResult.first == pointCloud1) {
            handleSelection(pointCloud1, selectedVertices1, vertexColors1);
        }
        else if (pickResult.first == pointCloud2) {
            handleSelection(pointCloud2, selectedVertices2, vertexColors2);
        }
    }

    // Update lines between matching selected vertices
    updateLineConnections(V1, V2_offset, selectedVertices1, selectedVertices2, radius_default / 3);

    std::cout << "Selected vertices in Point Cloud 1: ";
    for (int v : selectedVertices1) std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "Selected vertices in Point Cloud 2: ";
    for (int v : selectedVertices2) std::cout << v << " ";
    std::cout << std::endl;
}





// Rotation function
Eigen::MatrixXd rotatePointCloud(const Eigen::MatrixXd& V, double angle, const Eigen::Vector3d& axis) {
    // Calculate the centroid of the point cloud
    Eigen::Vector3d centroid = V.colwise().mean();

    // Create a rotation matrix
    Eigen::AngleAxisd rotation(angle, axis.normalized());  // 'angle' should be in radians
    Eigen::Matrix3d rotationMatrix = rotation.toRotationMatrix();

    // Shift the point cloud to the origin, rotate, and then shift back
    Eigen::MatrixXd V_rotated = V.rowwise() - centroid.transpose(); // Center at origin
    V_rotated = (rotationMatrix * V_rotated.transpose()).transpose(); // Rotate
    V_rotated.rowwise() += centroid.transpose(); // Shift back

    return V_rotated;
}


Eigen::VectorXd unitParameterizeBetweenPoints(const Eigen::MatrixXd& V, int startIdx, int endIdx) {
    // Ensure valid indices
    if (startIdx < 0 || startIdx >= V.rows() || endIdx < 0 || endIdx >= V.rows()) {
        throw std::invalid_argument("Invalid control point indices.");
    }

    bool isReversed = (startIdx > endIdx);
    int numVertices;
    Eigen::VectorXd params;

    if (isReversed) {
        // Wrap-around case: From startIdx -> end of V, and from start of V -> endIdx
        numVertices = V.rows() - startIdx + endIdx + 1;
        params.resize(numVertices);

        params(0) = 0.0; // First point is always 0
        params(numVertices - 1) = 1.0; // Last point is always 1

        double lengthCurve = 0.0;
        std::vector<double> segments(numVertices - 1, 0.0);

        // Calculate segment lengths for the wrap-around case
        // 1. From startIdx to end of V
        for (int i = 0; i < V.rows() - startIdx - 1; ++i) {
            segments[i] = distance2d(V.row(startIdx + i), V.row(startIdx + i + 1));
            lengthCurve += segments[i];
        }
        // 2. From the last point of V to the start of V
        segments[V.rows() - startIdx - 1] = distance2d(V.row(V.rows() - 1), V.row(0));
        lengthCurve += segments[V.rows() - startIdx - 1];

        // 3. From start of V to endIdx
        for (int i = 0; i < endIdx; ++i) {
            segments[V.rows() - startIdx + i] = distance2d(V.row(i), V.row(i + 1));
            lengthCurve += segments[V.rows() - startIdx + i];
        }

        // Compute parameter values
        for (int i = 1; i < numVertices - 1; ++i) {
            params(i) = params(i - 1) + (segments[i - 1] / lengthCurve);
        }
    }
    else {
        // Standard case: startIdx <= endIdx
        numVertices = endIdx - startIdx + 1;
        params.resize(numVertices);

        params(0) = 0.0; // First point is always 0
        params(numVertices - 1) = 1.0; // Last point is always 1

        double lengthCurve = 0.0;
        std::vector<double> segments(numVertices - 1, 0.0);

        // Calculate segment lengths between vertices
        for (int i = 0; i < numVertices - 1; ++i) {
            segments[i] = distance2d(V.row(startIdx + i), V.row(startIdx + i + 1));
            lengthCurve += segments[i];
        }

        // Compute parameter values for each vertex between startIdx and endIdx
        for (int i = 1; i < numVertices - 1; ++i) {
            params(i) = params(i - 1) + (segments[i - 1] / lengthCurve);
        }
    }

    return params;
}

Eigen::MatrixXd downsampleBetweenPoints(const Eigen::MatrixXd& V, int startIdx, int endIdx, int numSamples) {
    Eigen::MatrixXd downsampledPoints(numSamples, V.cols()); // Initialize a matrix to store downsampled points

    // Check if start and end indices are valid
    assert(startIdx >= 0 && startIdx < V.rows());
    assert(endIdx >= 0 && endIdx < V.rows());

    // Calculate total distance between the start and end points
    double totalDistance = 0.0;
    for (int i = startIdx; i != endIdx; i = (i + 1) % V.rows()) {
        totalDistance += (V.row((i + 1) % V.rows()) - V.row(i)).norm();
        if (i == endIdx) break;
    }

    // Calculate the distance increment for each sample point
    double distanceIncrement = totalDistance / (numSamples - 1);

    downsampledPoints.row(0) = V.row(startIdx); // Set the first sample point

    int currentIdx = startIdx;
    double accumulatedDistance = 0.0;

    for (int i = 1; i < numSamples; ++i) {
        double targetDistance = i * distanceIncrement;
        while (accumulatedDistance + (V.row((currentIdx + 1) % V.rows()) - V.row(currentIdx)).norm() < targetDistance) {
            accumulatedDistance += (V.row((currentIdx + 1) % V.rows()) - V.row(currentIdx)).norm();
            currentIdx = (currentIdx + 1) % V.rows();
            if (currentIdx == endIdx) break;
        }
        // Linear interpolation between current and next point
        double remainingDistance = targetDistance - accumulatedDistance;
        Eigen::VectorXd direction = (V.row((currentIdx + 1) % V.rows()) - V.row(currentIdx)).normalized();

        // Ensure dimension consistency
        assert(direction.size() == V.row(currentIdx).size());

        downsampledPoints.row(i) = V.row(currentIdx) + remainingDistance * direction;
    }

    return downsampledPoints;
}

Eigen::VectorXd adjustParamV2ToMatchV1(const Eigen::VectorXd& paramV1, const Eigen::VectorXd& paramV2, Eigen::VectorXd& alphas, Eigen::VectorXd& alphas_indices) {
    int sizeV1 = paramV1.size();
    int sizeV2 = paramV2.size();
    Eigen::VectorXd adjusted(sizeV1);

    // Ensure paramV2 has enough points for interpolation
    if (sizeV2 < 2) {
        throw std::invalid_argument("paramV2 must have at least 2 points for interpolation.");
    }

    for (int i = 0; i < sizeV1; ++i) {
        bool found = false; // Flag to check if we found the correct interval

        for (int j = 0; j < sizeV2 - 1; ++j) {
            // Check if paramV1[i] is within the interval (paramV2[j], paramV2[j+1])
            if (paramV2(j) <= paramV1(i) && paramV1(i) < paramV2(j + 1)) {
                // Calculate alpha
                alphas(i) = (paramV1(i) - paramV2(j)) / (paramV2(j + 1) - paramV2(j));
                alphas_indices(i) = j;
                // Interpolate adjusted value using the calculated alpha
                adjusted(i) = (1 - alphas(i)) * paramV2(j) + alphas(i) * paramV2(j + 1);

                found = true; // Mark that we found the interval
                break; // Exit the inner loop since we found the interval
            }
        }

        if (!found) {
            // If no valid interval was found, handle accordingly (e.g., set adjusted to a default value or throw an error)
            adjusted(i) = paramV2(0); // Default to the first element, or you could throw an error
            alphas(i) = 0.0; // Default alpha
            alphas_indices(i) = 0.0; // Default alpha
        }
    }

    return adjusted;
}






void parameterizeWithControls(const Eigen::MatrixXd& V1, const Eigen::MatrixXd& V2,
    const std::vector<int>& selectedVertices1,
    const std::vector<int>& selectedVertices2) {
    std::cout << "Parameterization with landmarks called." << std::endl;

    int numLandmarks1 = selectedVertices1.size();
    int numLandmarks2 = selectedVertices2.size();
    int numLandmarks = std::min(numLandmarks1, numLandmarks2);

    if (numLandmarks < 2) {
        std::cerr << "Not enough landmarks to parameterize. At least 2 landmarks are required on both curves." << std::endl;
        return; // Exit if there are not enough landmarks
    }

    // For each valid pair of consecutive landmarks, compute the parameterization for both curves
    for (int i = 0; i < numLandmarks; ++i) {
        int startIdx1 = selectedVertices1[i];
        int endIdx1 = (i + 1 == numLandmarks) ? selectedVertices1[0] : selectedVertices1[i + 1];

        int startIdx2 = selectedVertices2[i];
        int endIdx2 = (i + 1 == numLandmarks) ? selectedVertices2[0] : selectedVertices2[i + 1];

        // Ensure the indices are valid
        if (startIdx1 >= V1.rows() || endIdx1 >= V1.rows() || startIdx2 >= V2.rows() || endIdx2 >= V2.rows()) {
            std::cerr << "Invalid indices for parameterization." << std::endl;
            return;
        }

        // Get sub-vectors for V1 and V2 based on the selected landmarks
        Eigen::MatrixXd V1_sub, V2_sub;

        // For V1
        if (endIdx1 >= startIdx1) {
            // Normal case
            V1_sub = V1.block(startIdx1, 0, endIdx1 - startIdx1 + 1, V1.cols());
        }
        else {
            // Wrap-around case
            V1_sub = V1.block(startIdx1, 0, V1.rows() - startIdx1, V1.cols()); // From startIdx1 to end
            Eigen::MatrixXd wrapBlock = V1.block(0, 0, endIdx1 + 1, V1.cols()); // From start to endIdx1
            V1_sub.conservativeResize(V1_sub.rows() + wrapBlock.rows(), V1.cols());
            V1_sub.bottomRows(wrapBlock.rows()) = wrapBlock; // Append wrapBlock
        }

        // For V2
        if (endIdx2 >= startIdx2) {
            // Normal case
            V2_sub = V2.block(startIdx2, 0, endIdx2 - startIdx2 + 1, V2.cols());
        }
        else {
            // Wrap-around case
            V2_sub = V2.block(startIdx2, 0, V2.rows() - startIdx2, V2.cols()); // From startIdx2 to end
            Eigen::MatrixXd wrapBlock = V2.block(0, 0, endIdx2 + 1, V2.cols()); // From start to endIdx2
            V2_sub.conservativeResize(V2_sub.rows() + wrapBlock.rows(), V2.cols());
            V2_sub.bottomRows(wrapBlock.rows()) = wrapBlock; // Append wrapBlock
        }

        //// Get sub-vectors for V1 and V2 based on the selected landmarks
        //Eigen::MatrixXd V1_sub = V1.block(startIdx1, 0, endIdx1 - startIdx1 + 1, V1.cols());
        //Eigen::MatrixXd V2_sub = V2.block(startIdx2, 0, endIdx2 - startIdx2 + 1, V2.cols());

        // Get equalized point clouds
        auto [V1_equalized, V2_equalized] = getEqualizedPointClouds(V1_sub, V2_sub);

        // Print V1_equalized and V2_equalized vertices
        std::cout << "V1_equalized vertices:\n" << V1_equalized << std::endl;
        std::cout << "V2_equalized vertices:\n" << V2_equalized << std::endl;

        // Unit parameterization for equalized point clouds
        Eigen::VectorXd paramV1 = unitParameterizeBetweenPoints(V1_equalized, 0, V1_equalized.rows() - 1);
        Eigen::VectorXd paramV2 = unitParameterizeBetweenPoints(V2_equalized, 0, V2_equalized.rows() - 1);

        Eigen::VectorXd alphas(paramV1.size());
        Eigen::VectorXd alphas_indices(paramV1.size());
        // Interpolate paramV2 to match the size of paramV1
        Eigen::VectorXd interpolatedV2 = adjustParamV2ToMatchV1(paramV1, paramV2, alphas, alphas_indices);

        // Output the results
        std::cout << "Parameterization for V1:\n" << paramV1 << std::endl;
        std::cout << "Interpolated V2:\n" << interpolatedV2 << std::endl;
        std::cout << "alphas V2:\n" << alphas << std::endl;
        std::cout << "alphas V2 indices:\n" << alphas_indices << std::endl;

        // Create NewV2 by interpolating between vertices of V2_equalized using alphas and indices
        Eigen::MatrixXd NewV2(interpolatedV2.size(), V2_equalized.cols());
        for (int j = 0; j < interpolatedV2.size(); ++j) {
            int idx = static_cast<int>(alphas_indices[j]);  // Get the base index
            double alpha = alphas[j];  // Get the interpolation factor

            // Ensure idx is valid for interpolation
            if (idx + 1 < V2_equalized.rows()) {
                // Linear interpolation between V2_equalized[idx] and V2_equalized[idx + 1]
                NewV2.row(j) = (1.0 - alpha) * V2_equalized.row(idx) + alpha * V2_equalized.row(idx + 1);
            }
            else {
                // Handle boundary condition where idx + 1 would be out of bounds
                NewV2.row(j) = V2_equalized.row(idx);
            }
        }


        // Register point clouds
        polyscope::registerPointCloud("V1 Equalized___" + std::to_string(i), V1_equalized);
        //polyscope::registerPointCloud("New V2", NewV2);
        //polyscope::registerPointCloud("V2 Equalized", V2_equalized);
        //polyscope::registerPointCloud("V1 Sub", V1_sub);
        //polyscope::registerPointCloud("V2 Sub", V2_sub);
        polyscope::registerPointCloud("V2 new___" + std::to_string(i), NewV2);


        // Output NewV2 for inspection
        //std::cout << "NewV2 vertices:\n" << NewV2 << std::endl;

        // Break after the first pair of landmarks (to be adjusted as needed)
        //break;
    }

    std::cout << "Landmark-based parameterization complete (with wrap-around handling)." << std::endl;
}



// Main function with slider for rotation of both point clouds
void showSideBySideSelectionWithVertexSelection(const Eigen::MatrixXd& V1, const Eigen::MatrixXd& V2) {
    polyscope::init();

    double radius_default = 0.0005 * 3;

    std::vector<int> selectedVertices1, selectedVertices2;
    std::vector<std::array<double, 3>> vertexColors1, vertexColors2;

    // Initial rotation angles (in radians) for both point clouds
    static double rotationAngleRadians1 = 0.0;
    static double rotationAngleRadians2 = 0.0;

    // Axis of rotation (around Z-axis)
    Eigen::Vector3d rotationAxis(0, 0, 1);

    // Register the first point cloud
    auto* pointCloud1 = registerPointCloudWithColors("Point Cloud 1", V1, radius_default, vertexColors1);

    // Calculate and adjust offsets for V2
    Eigen::MatrixXd V2_offset = calculateAndAdjustOffsets(V1, V2);

    // Register the second point cloud (rotated)
    auto* pointCloud2 = registerPointCloudWithColors("Point Cloud 2", rotatePointCloud(V2_offset, rotationAngleRadians2, rotationAxis), radius_default, vertexColors2);

    // Set user callback for handling UI and updates
    polyscope::state::userCallback = [&]() {
        // Slider for rotating the first model (V1)
        float rotationAngleDegrees1 = static_cast<float>(rotationAngleRadians1 * 180.0 / M_PI);  // Convert radians to degrees
        if (ImGui::SliderFloat("Rotation Angle Left Model (Degrees)", &rotationAngleDegrees1, -180.0f, 180.0f)) {
            // Convert degrees back to radians
            rotationAngleRadians1 = rotationAngleDegrees1 * M_PI / 180.0;
        }

        // Slider for rotating the second model (V2)
        float rotationAngleDegrees2 = static_cast<float>(rotationAngleRadians2 * 180.0 / M_PI);  // Convert radians to degrees
        if (ImGui::SliderFloat("Rotation Angle Right Model (Degrees)", &rotationAngleDegrees2, -180.0f, 180.0f)) {
            // Convert degrees back to radians
            rotationAngleRadians2 = rotationAngleDegrees2 * M_PI / 180.0;
        }

        // Update the point clouds with the new rotation angles
        Eigen::MatrixXd V1_rotated = rotatePointCloud(V1, rotationAngleRadians1, rotationAxis); // Rotate the first point cloud
        pointCloud1->updatePointPositions(V1_rotated); // Update Polyscope visualization for the first model

        Eigen::MatrixXd V2_rotated = rotatePointCloud(V2_offset, rotationAngleRadians2, rotationAxis); // Rotate the second point cloud
        pointCloud2->updatePointPositions(V2_rotated); // Update Polyscope visualization for the second model

        // Handle user selection of vertices
        handleUserSelection(pointCloud1, pointCloud2, V1_rotated, V2_rotated, selectedVertices1, selectedVertices2, vertexColors1, vertexColors2, radius_default);

        // Add a button for "parameterize"
        if (ImGui::Button("parameterize")) {
            // Call the parameterize function when the button is pressed
            parameterizeWithControls(V1, V2_offset, selectedVertices1, selectedVertices2);
        }

        // Add a button to reset selected points
        if (ImGui::Button("Reset Selections")) {
            // Clear selections
            selectedVertices1.clear();
            selectedVertices2.clear();
            vertexColors1.clear();
            vertexColors2.clear();

            // Reinitialize vertex colors to the default state
            vertexColors1.resize(V1_rotated.rows(), { 1.0, 1.0, 1.0 }); // Set default color (white) for each vertex in V1
            vertexColors2.resize(V2_rotated.rows(), { 1.0, 1.0, 1.0 }); // Set default color (white) for each vertex in V2

            // Update the point cloud colors to reflect no selection
            pointCloud1->updatePointPositions(V1_rotated);  // Reset colors for V1
            pointCloud2->updatePointPositions(V2_rotated);  // Reset colors for V2

            // Remove the existing curve network to reset the lines
            if (polyscope::hasCurveNetwork("Corresponding Lines")) {
                polyscope::removeCurveNetwork("Corresponding Lines"); // Reset the existing curve network
            }

            std::cout << "Selections reset." << std::endl;
        }

    };

    polyscope::show();
}



