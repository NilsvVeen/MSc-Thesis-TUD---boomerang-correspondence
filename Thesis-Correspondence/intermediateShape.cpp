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
#include <polyscope/surface_mesh.h>

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

// Main
void main_phase2(Eigen::MatrixXd V1, Eigen::MatrixXi F1, Eigen::MatrixXd V2, Eigen::MatrixXi F2 ) {
    Eigen::MatrixXd V_new;


    // Interpolate vertices
    double t = 0.5; // Intermediate shape
    V_new = interpolateVertices(V1, V2, t);

    // Optimize for rigidity
    V_new = performARAP(V_new, F1);

    // Visualize with Polyscope
    polyscope::init();
    polyscope::registerSurfaceMesh("Interpolated Shape", V_new, F1);
    polyscope::registerSurfaceMesh("Shape 1", V1, F1);
    polyscope::registerSurfaceMesh("Shape 2", V2, F2);
    polyscope::show();

}
