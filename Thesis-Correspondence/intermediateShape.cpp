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

// Global variables for mesh storage
Eigen::MatrixXd V1, V2, V_new;
Eigen::MatrixXi F1, F2;
double t = 0.5; // Interpolation factor
Eigen::Matrix3d rotationMatrix1 = Eigen::Matrix3d::Identity();
Eigen::Matrix3d rotationMatrix2 = Eigen::Matrix3d::Identity();
Eigen::Vector2d translation1(0, 0);
Eigen::Vector2d translation2(0, 0);

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

// Interpolate vertices
Eigen::MatrixXd interpolateVertices(const Eigen::MatrixXd& V1, const Eigen::MatrixXd& V2, double t) {
    return (1.0 - t) * V1 + t * V2;
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

    // Interpolate vertices
    V_new = interpolateVertices(V1, V2, t);

    // Optimize for rigidity
    V_new = performARAP(V_new, F1);

    // Update Polyscope visualization
    polyscope::getSurfaceMesh("Interpolated Shape")->updateVertexPositions(V_new);
}

// Update transformations and Polyscope
void updateMeshTransformations() {
    Eigen::MatrixXd transformedV1 = applyTransformations(V1, rotationMatrix1, translation1);
    Eigen::MatrixXd transformedV2 = applyTransformations(V2, rotationMatrix2, translation2);

    polyscope::getSurfaceMesh("Shape 1")->updateVertexPositions(transformedV1);
    polyscope::getSurfaceMesh("Shape 2")->updateVertexPositions(transformedV2);
}

// Main function
void main_phase2(Eigen::MatrixXd inputV1, Eigen::MatrixXi inputF1, Eigen::MatrixXd inputV2, Eigen::MatrixXi inputF2) {
    V1 = inputV1;
    F1 = inputF1;
    V2 = inputV2;
    F2 = inputF2;

    polyscope::init();

    polyscope::registerSurfaceMesh("Shape 1", V1, F1);
    polyscope::registerSurfaceMesh("Shape 2", V2, F2);

    V_new = interpolateVertices(V1, V2, t);
    polyscope::registerSurfaceMesh("Interpolated Shape", V_new, F1);

    polyscope::state::userCallback = []() {
        ImGui::Text("Object 1 Transformations");
        static float rotation1 = 0.0f;
        static float translationX1 = 0.0f;
        static float translationY1 = 0.0f;

        if (ImGui::SliderFloat("Rotation 1 (degrees)", &rotation1, 0.0f, 360.0f)) {
            double angleRad = rotation1 * M_PI / 180.0;
            rotationMatrix1 << cos(angleRad), -sin(angleRad), 0,
                sin(angleRad), cos(angleRad), 0,
                0, 0, 1;
            updateMeshTransformations();
        }
        if (ImGui::SliderFloat("Translation X1", &translationX1, -50.0f, 50.0f)) {
            translation1.x() = translationX1;
            updateMeshTransformations();
        }
        if (ImGui::SliderFloat("Translation Y1", &translationY1, -50.0f, 50.0f)) {
            translation1.y() = translationY1;
            updateMeshTransformations();
        }

        ImGui::Separator();

        ImGui::Text("Object 2 Transformations");
        static float rotation2 = 0.0f;
        static float translationX2 = 0.0f;
        static float translationY2 = 0.0f;

        if (ImGui::SliderFloat("Rotation 2 (degrees)", &rotation2, 0.0f, 360.0f)) {
            double angleRad = rotation2 * M_PI / 180.0;
            rotationMatrix2 << cos(angleRad), -sin(angleRad), 0,
                sin(angleRad), cos(angleRad), 0,
                0, 0, 1;
            updateMeshTransformations();
        }
        if (ImGui::SliderFloat("Translation X2", &translationX2, -50.0f, 50.0f)) {
            translation2.x() = translationX2;
            updateMeshTransformations();
        }
        if (ImGui::SliderFloat("Translation Y2", &translationY2, -50.0f, 50.0f)) {
            translation2.y() = translationY2;
            updateMeshTransformations();
        }

        ImGui::Separator();

        ImGui::Text("Interpolation");
        float t_float = static_cast<float>(t);
        if (ImGui::SliderFloat("Interpolation Factor (t)", &t_float, 0.0f, 1.0f)) {
            t = static_cast<double>(t_float);
        }
        if (ImGui::Button("Generate New Shape")) {
            updateIntermediateShape();
        }
    };

    polyscope::show();
}
