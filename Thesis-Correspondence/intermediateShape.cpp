// Undefine PI before any other includes
#ifdef PI
#undef PI
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

// Main function
void main_phase2(Eigen::MatrixXd inputV1, Eigen::MatrixXi inputF1, Eigen::MatrixXd inputV2, Eigen::MatrixXi inputF2) {
    // Initialize global variables
    V1 = inputV1;
    F1 = inputF1;
    V2 = inputV2;
    F2 = inputF2;

    // Initialize Polyscope
    polyscope::init();

    // Register initial meshes
    polyscope::registerSurfaceMesh("Shape 1", V1, F1);
    polyscope::registerSurfaceMesh("Shape 2", V2, F2);

    // Register the interpolated shape mesh
    V_new = interpolateVertices(V1, V2, t); // Initial interpolation
    polyscope::registerSurfaceMesh("Interpolated Shape", V_new, F1);

    // Add GUI callback to update the intermediate shape
    polyscope::state::userCallback = []() {
        ImGui::Text("Update Intermediate Shape");

        // Use float for the slider since ImGui does not support double
        float t_float = static_cast<float>(t);
        if (ImGui::SliderFloat("Interpolation Factor (t)", &t_float, 0.0f, 1.0f)) {
            t = static_cast<double>(t_float); // Update the global double variable
        }
        if (ImGui::Button("Generate New Shape")) {
            updateIntermediateShape();
        }
    };

    // Show the Polyscope GUI
    polyscope::show();
}

