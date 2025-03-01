// Undefine PI before any other includes
#ifdef PI
#undef PI
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Libigl includes
#include <igl/PI.h>
#include <igl/arap.h>

// Other includes
#include <Eigen/Dense>
#include <Polyscope/polyscope.h>
#include <Polyscope/surface_mesh.h>
#include <iostream>
#include <vector>
#include "stl_utils.h"
#include "file_utils.h"
#include "file_utils.h"

// Global variables for mesh storage
std::vector<Eigen::MatrixXd> vertices;
std::vector<Eigen::MatrixXi> faces;
std::vector<Eigen::Matrix3d> rotationMatrices;
std::vector<Eigen::Vector3d> translations;
Eigen::MatrixXd V_new;
Eigen::MatrixXi F_new;
std::vector<double> weights; // Interpolation weights
std::string outputFolderCopy;

// Function to compute the center of gravity of a mesh
Eigen::RowVector3d computeCOG(const Eigen::MatrixXd& V) {
    return V.colwise().mean();
}

// Apply transformations (rotation and translation) around COG
Eigen::MatrixXd applyTransformations(const Eigen::MatrixXd& V, const Eigen::Matrix3d& R, const Eigen::Vector3d& T) {
    Eigen::RowVector3d cog = computeCOG(V);
    Eigen::MatrixXd transformedV = V;

    // Translate to origin (COG), apply rotation, and translate back
    transformedV.rowwise() -= cog;
    transformedV.leftCols<2>() = (transformedV.leftCols<2>() * R.topLeftCorner<2, 2>().transpose());
    transformedV.rowwise() += cog;

    // Apply translation in XY plane
    transformedV.col(0).array() += T(0);
    transformedV.col(1).array() += T(1);
    transformedV.col(2).array() += T(2);

    return transformedV;
}

// Interpolate vertices across multiple shapes
Eigen::MatrixXd interpolateVertices(const std::vector<Eigen::MatrixXd>& vertices, const std::vector<double>& weights) {
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(vertices[0].rows(), vertices[0].cols());

    for (size_t i = 0; i < vertices.size(); ++i) {
        result += weights[i] * vertices[i];
    }

    return result;
}

// Perform ARAP for rigidity preservation
Eigen::MatrixXd performARAP(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    igl::ARAPData arap_data;
    arap_data.with_dynamics = false;

    Eigen::VectorXi b(1); // Fix one vertex
    b(0) = 0;
    Eigen::MatrixXd bc(1, 3);
    bc = V.row(0);

    igl::arap_precomputation(V, F, 3, b, arap_data);
    Eigen::MatrixXd V_new = V;
    igl::arap_solve(bc, arap_data, V_new);

    return V_new;
}

// Callback to update the intermediate shape
void updateIntermediateShape() {
    std::cout << "Updating intermediate shape..." << std::endl;

    std::vector<Eigen::MatrixXd> transformedVertices;
    for (size_t i = 0; i < vertices.size(); ++i) {
        transformedVertices.push_back(applyTransformations(vertices[i], rotationMatrices[i], translations[i]));
    }

    // Interpolate vertices
    V_new = interpolateVertices(transformedVertices, weights);

    // Optimize for rigidity
    V_new = performARAP(V_new, F_new);

    // Update Polyscope visualization
    polyscope::getSurfaceMesh("Interpolated Shape")->updateVertexPositions(V_new);
}

// Update transformations and Polyscope
void updateMeshTransformations() {
    for (size_t i = 0; i < vertices.size(); ++i) {
        Eigen::MatrixXd transformedV = applyTransformations(vertices[i], rotationMatrices[i], translations[i]);
        polyscope::getSurfaceMesh("Shape " + std::to_string(i + 1))->updateVertexPositions(transformedV);
    }
}

Eigen::VectorXd selectSubsetIndices(int totalVertices, int numSelected) {
    Eigen::VectorXd indices = Eigen::VectorXd::Constant(totalVertices, -1); // Default to -1 (not selected)

    for (int i = 0; i < numSelected; ++i) {
        int idx = i * (totalVertices / numSelected); // Evenly spaced selection
        indices(idx) = double(i) / (numSelected - 1); // Normalize to [0,1]
    }

    return indices;
}


#include <Eigen/Core>
#include <queue>
#include <vector>
#include <limits>

// Function to compute geodesic distances from a starting vertex
Eigen::VectorXd computeSortedVertexIndices(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int startIdx = 0) {
    int numVertices = V.rows();

    // Priority queue for Dijkstra's algorithm (min-heap based on distance)
    using Pair = std::pair<double, int>;  // (distance, vertex index)
    std::priority_queue<Pair, std::vector<Pair>, std::greater<Pair>> pq;

    // Distance array initialized to infinity
    std::vector<double> distances(numVertices, std::numeric_limits<double>::infinity());
    distances[startIdx] = 0.0;

    // Result vector (stores indices sorted by increasing distance)
    std::vector<int> sortedIndices;

    // Push the starting vertex into the queue
    pq.emplace(0.0, startIdx);

    // Dijkstra’s algorithm for shortest path in an unstructured mesh
    while (!pq.empty()) {
        double currDist = pq.top().first;
        int currIdx = pq.top().second;
        pq.pop();

        // If already visited, skip
        if (std::find(sortedIndices.begin(), sortedIndices.end(), currIdx) != sortedIndices.end()) continue;

        // Add to sorted indices
        sortedIndices.push_back(currIdx);

        // Iterate over neighboring vertices (connected via faces)
        for (int i = 0; i < F.rows(); ++i) {
            for (int j = 0; j < 3; ++j) {
                if (F(i, j) == currIdx) {
                    // Neighboring vertex
                    int neighborIdx = F(i, (j + 1) % 3);

                    // Compute Euclidean distance
                    double newDist = currDist + (V.row(currIdx) - V.row(neighborIdx)).norm();

                    // Relaxation step
                    if (newDist < distances[neighborIdx]) {
                        distances[neighborIdx] = newDist;
                        pq.emplace(newDist, neighborIdx);
                    }
                }
            }
        }
    }

    // Convert to Eigen::VectorXd
    Eigen::VectorXd sortedEigenIndices(sortedIndices.size());
    for (int i = 0; i < sortedIndices.size(); ++i) {
        sortedEigenIndices(i) = sortedIndices[i];
    }

    return sortedEigenIndices;
}


// Compute shortest path distances over the mesh using Dijkstra's algorithm
Eigen::VectorXd computeGeodesicDistances(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int startIdx = 0) {
    int numVertices = V.rows();

    // Priority queue for Dijkstra (min-heap based on distance)
    using Pair = std::pair<double, int>;  // (distance, vertex index)
    std::priority_queue<Pair, std::vector<Pair>, std::greater<Pair>> pq;

    // Distance array initialized to infinity
    std::vector<double> distances(numVertices, std::numeric_limits<double>::infinity());
    distances[startIdx] = 0.0;  // Distance to itself is 0

    // Push the starting vertex into the queue
    pq.emplace(0.0, startIdx);

    // Adjacency list for edge-based traversal
    std::vector<std::vector<std::pair<int, double>>> adjacency(numVertices);
    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int v1 = F(i, j);
            int v2 = F(i, (j + 1) % 3);
            int v3 = F(i, (j + 2) % 3);

            double d12 = (V.row(v1) - V.row(v2)).norm();
            double d13 = (V.row(v1) - V.row(v3)).norm();
            double d23 = (V.row(v2) - V.row(v3)).norm();

            adjacency[v1].emplace_back(v2, d12);
            adjacency[v1].emplace_back(v3, d13);
            adjacency[v2].emplace_back(v1, d12);
            adjacency[v2].emplace_back(v3, d23);
            adjacency[v3].emplace_back(v1, d13);
            adjacency[v3].emplace_back(v2, d23);
        }
    }

    // Dijkstra's algorithm for geodesic distances
    while (!pq.empty()) {
        double currDist = pq.top().first;
        int currIdx = pq.top().second;
        pq.pop();

        // Iterate over neighboring vertices
        for (const auto& neighbor : adjacency[currIdx]) {
            int neighborIdx = neighbor.first;
            double edgeLength = neighbor.second;
            double newDist = currDist + edgeLength;

            // Relaxation step
            if (newDist < distances[neighborIdx]) {
                distances[neighborIdx] = newDist;
                pq.emplace(newDist, neighborIdx);
            }
        }
    }

    // Convert distances to Eigen::VectorXd
    Eigen::VectorXd distanceVector(numVertices);
    for (int i = 0; i < numVertices; ++i) {
        distanceVector(i) = distances[i];
    }

    return distanceVector;
}


// Main function
void main_phase2(const std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>>& inputShapes, const std::string outputFolder
) {

    outputFolderCopy = outputFolder;
    vertices.clear();
    faces.clear();
    rotationMatrices.clear();
    translations.clear();
    weights.clear();

    for (const auto& shape : inputShapes) {
        vertices.push_back(shape.first);
        faces.push_back(shape.second);
        rotationMatrices.emplace_back(Eigen::Matrix3d::Identity());
        translations.emplace_back(Eigen::Vector3d::Zero());
        weights.push_back(1.0 / inputShapes.size()); // Default equal weights
    }

    F_new = faces[0]; // Assume all shapes have the same topology

    polyscope::init();

    //int numSelected = 5000; // Only highlight 100 points
    //Eigen::VectorXd indices = selectSubsetIndices(vertices[0].rows(), numSelected);
    //Eigen::VectorXd indices = computeSortedVertexIndices(vertices[0], F_new);
    Eigen::VectorXd indices = computeGeodesicDistances(vertices[0], F_new);

    std::cout << "indices type 3" << std::endl;

    for (size_t i = 0; i < vertices.size(); ++i) {
        //auto* obj = polyscope::registerPointCloud("Shape " + std::to_string(i + 1) + "-Points", vertices[i]);
        auto* obj2 = polyscope::registerSurfaceMesh("Shape " + std::to_string(i + 1), vertices[i], faces[i]);
        //obj->addScalarQuantity("color", indices);
        obj2->addVertexScalarQuantity("color", indices);
    }

    // Interpolate the transformed meshes
    V_new = interpolateVertices(vertices, weights);
    auto* objInterpolated = polyscope::registerSurfaceMesh("Interpolated Shape", V_new, F_new);
    objInterpolated->addVertexScalarQuantity("color", indices); // Add color to interpolated shape

    polyscope::state::userCallback = []() {
        for (size_t i = 0; i < vertices.size(); ++i) {
            ImGui::Text("Object %zu Transformations", i + 1);
            static float rotation = 0.0f;
            static float translationX = 0.0f;
            static float translationY = 0.0f;

            if (ImGui::SliderFloat(("Rotation " + std::to_string(i + 1) + " (degrees)").c_str(), &rotation, 0.0f, 360.0f)) {
                double angleRad = rotation * M_PI / 180.0;
                rotationMatrices[i] << cos(angleRad), -sin(angleRad), 0,
                    sin(angleRad), cos(angleRad), 0,
                    0, 0, 1;
                updateMeshTransformations();
            }
            if (ImGui::SliderFloat(("Translation X" + std::to_string(i + 1)).c_str(), &translationX, -500.0f, 500.0f)) {
                translations[i].x() = translationX;
                updateMeshTransformations();
            }
            if (ImGui::SliderFloat(("Translation Y" + std::to_string(i + 1)).c_str(), &translationY, -500.0f, 500.0f)) {
                translations[i].y() = translationY;
                updateMeshTransformations();
            }

            ImGui::Separator();
        }

        ImGui::Text("Interpolation");
        static std::vector<float> weightsFloat;
        if (weightsFloat.empty()) {
            weightsFloat.resize(weights.size(), 1.0f / weights.size());
        }

        // Track the sum of weights
        float totalWeight = 0.0f;
        for (float w : weightsFloat) {
            totalWeight += w;
        }

        for (size_t i = 0; i < weights.size(); ++i) {
            float oldValue = weightsFloat[i];
            if (ImGui::SliderFloat(("Weight " + std::to_string(i + 1)).c_str(), &weightsFloat[i], 0.0f, 1.0f)) {
                // Adjust other weights to keep total <= 1
                float diff = weightsFloat[i] - oldValue;
                float remainingAdjustment = diff;

                // Redistribute the adjustment proportionally across other sliders
                for (size_t j = 0; j < weightsFloat.size(); ++j) {
                    if (j == i) continue; // Skip the currently modified slider
                    float adjustment = -remainingAdjustment * (weightsFloat[j] / (totalWeight - oldValue));
                    weightsFloat[j] = std::max(0.0f, std::min(1.0f, weightsFloat[j] + adjustment));
                    remainingAdjustment += adjustment; // Keep track of remaining adjustment
                }

                // Recalculate total weight and ensure normalization
                totalWeight = 0.0f;
                for (float w : weightsFloat) {
                    totalWeight += w;
                }

                // Normalize all sliders if total exceeds 1
                if (totalWeight > 1.0f) {
                    for (float& w : weightsFloat) {
                        w /= totalWeight;
                    }
                }

                // Update the double weights
                for (size_t k = 0; k < weights.size(); ++k) {
                    weights[k] = static_cast<double>(weightsFloat[k]);
                }
            }
        }

        //// Button to calculate the average position of the first model and translate all other models to align with it
        if (ImGui::Button("Align All Models to First Model's Average Position")) {
            // Calculate the average position of the first model
            std::cout << "1" << std::endl;
            Eigen::Vector3d firstModelAverage = vertices[0].colwise().mean().transpose();

            std::cout << "2222" << std::endl;


            // Translate all other models to align with the first model's average position
            for (size_t i = 1; i < vertices.size(); ++i) { // Skip the first model

                std::cout << "i" << std::endl;

                Eigen::Vector3d currentAverage = vertices[i].colwise().mean().transpose();
                Eigen::Vector3d translationVector = firstModelAverage - currentAverage;

                // Apply translation to the transformation vector
                translations[i].x() = translationVector.x();
                translations[i].y() = translationVector.y();
                translations[i].z() = translationVector.z();

                // Update mesh transformations
                updateMeshTransformations();
            }
        }



        //std::string outputFolder = "Matched/";
        //// Button to save all transformed shapes and the interpolated shape
        if (ImGui::Button("Save All Shapes")) {
            // Explicitly capture 'outputFolder' in the lambda function
                // Ensure the output folder exists and is clean
                createDirectory(outputFolderCopy);
                clearDirectory(outputFolderCopy);

                // Save each transformed shape
                for (size_t i = 0; i < vertices.size(); ++i) {
                    // Apply transformations to the vertices
                    Eigen::MatrixXd transformedVertices = (rotationMatrices[i] * vertices[i].transpose()).transpose();
                    transformedVertices.rowwise() += translations[i].transpose();

                    // Construct filename for each shape
                    std::string shapeFilename = outputFolderCopy + "/Shape_" + std::to_string(i + 1) + ".obj";

                    // Save the transformed vertices and faces
                    saveMeshToFile(shapeFilename, transformedVertices, faces[i]);
                }

                // Save the interpolated shape
                std::string interpolatedFilename = outputFolderCopy + "/Interpolated_Shape.obj";
                saveMeshToFile(interpolatedFilename, V_new, F_new);

                // Notify the user
                std::cout << "All shapes saved to folder: " << outputFolderCopy << std::endl;
        }







        ImGui::Separator();

        if (ImGui::Button("Generate New Shape")) {
            updateIntermediateShape();
        }
    };

    polyscope::show();
}

