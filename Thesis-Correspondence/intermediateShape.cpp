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

// Global variables for mesh storage
std::vector<Eigen::MatrixXd> vertices;
std::vector<Eigen::MatrixXi> faces;
std::vector<Eigen::Matrix3d> rotationMatrices;
std::vector<Eigen::Vector2d> translations;
Eigen::MatrixXd V_new;
Eigen::MatrixXi F_new;
std::vector<double> weights; // Interpolation weights

// Function to compute the center of gravity of a mesh
Eigen::RowVector3d computeCOG(const Eigen::MatrixXd& V) {
    return V.colwise().mean();
}

// Apply transformations (rotation and translation) around COG
Eigen::MatrixXd applyTransformations(const Eigen::MatrixXd& V, const Eigen::Matrix3d& R, const Eigen::Vector2d& T) {
    Eigen::RowVector3d cog = computeCOG(V);
    Eigen::MatrixXd transformedV = V;

    // Translate to origin (COG), apply rotation, and translate back
    transformedV.rowwise() -= cog;
    transformedV.leftCols<2>() = (transformedV.leftCols<2>() * R.topLeftCorner<2, 2>().transpose());
    transformedV.rowwise() += cog;

    // Apply translation in XY plane
    transformedV.col(0).array() += T(0);
    transformedV.col(1).array() += T(1);

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

// Main function
void main_phase2(const std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>>& inputShapes) {
    vertices.clear();
    faces.clear();
    rotationMatrices.clear();
    translations.clear();
    weights.clear();

    for (const auto& shape : inputShapes) {
        vertices.push_back(shape.first);
        faces.push_back(shape.second);
        rotationMatrices.emplace_back(Eigen::Matrix3d::Identity());
        translations.emplace_back(Eigen::Vector2d::Zero());
        weights.push_back(1.0 / inputShapes.size()); // Default equal weights
    }

    F_new = faces[0]; // Assume all shapes have the same topology

    polyscope::init();

    for (size_t i = 0; i < vertices.size(); ++i) {
        polyscope::registerSurfaceMesh("Shape " + std::to_string(i + 1), vertices[i], faces[i]);
    }

    // Interpolate the transformed meshes
    V_new = interpolateVertices(vertices, weights);
    polyscope::registerSurfaceMesh("Interpolated Shape", V_new, F_new);

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


        if (ImGui::Button("Generate New Shape")) {
            updateIntermediateShape();
        }
    };

    polyscope::show();
}
