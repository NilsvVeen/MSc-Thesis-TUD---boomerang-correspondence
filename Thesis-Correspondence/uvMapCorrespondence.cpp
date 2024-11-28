#include <Eigen/Core>       // For Eigen::MatrixXd, Eigen::MatrixXi, and basic matrix operations
#include <Eigen/Dense>      // For additional dense matrix operations (e.g., decomposition)
#include <Eigen/Geometry>   // For geometric transformations if required (e.g., rotations)
#include <iostream>
#include <Eigen/Dense>
#include <vector>

#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>
#include <polyscope/surface_mesh.h>





// Function to check if a point is inside a triangle in UV space
bool isPointInTriangle(const Eigen::RowVector2d& point, const Eigen::RowVector2d& A,
    const Eigen::RowVector2d& B, const Eigen::RowVector2d& C) {
    Eigen::RowVector2d v0 = C - A;
    Eigen::RowVector2d v1 = B - A;
    Eigen::RowVector2d v2 = point - A;

    double dot00 = v0.dot(v0);
    double dot01 = v0.dot(v1);
    double dot02 = v0.dot(v2);
    double dot11 = v1.dot(v1);
    double dot12 = v1.dot(v2);

    double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    return (u >= 0) && (v >= 0) && (u + v <= 1);
}

void UVToCorrespondence(
    const Eigen::MatrixXd& V1,  // Mesh 1 vertices       n x 3
    const Eigen::MatrixXi& F1, // Mesh 1 faces          m x 3
    const Eigen::MatrixXd& B1, // Boundary 1 vertices   p x 3
    const Eigen::MatrixXd& UV1, // UV map of mesh 1     n x 2

    const Eigen::MatrixXd& V2,  // Mesh 2 vertices       s x 3
    const Eigen::MatrixXi& F2,  // Mesh 2 faces          q x 3
    const Eigen::MatrixXd& B2,  // Boundary 2 vertices   p x 3
    const Eigen::MatrixXd& UV2  // UV map of mesh 2     s x 2
) {
    // Function implementation here.


    polyscope::init();
    polyscope::registerSurfaceMesh("M1", V1, F1);
    polyscope::registerSurfaceMesh("M2", V2, F2);

   // Point clouds for registration
    std::vector<Eigen::RowVector3d> pointCloudA; // Holds 3D points from V1
    std::vector<Eigen::RowVector3d> pointCloudB; // Holds 3D points from face vertices in V2

    for (int i = 0; i < V1.rows(); ++i) {
        Eigen::RowVector3d a = V1.row(i);     // Vertex in 3D
        Eigen::RowVector2d b = UV1.row(i);   // Corresponding UV coordinate

        // Find the face in UV2 that contains b
        //for (int j = 0; j < F2.rows(); ++j) {
        for (int j = F2.rows()-1; j > 1; --j) {
            // Get the UV coordinates of the face vertices in UV2
            Eigen::RowVector2d A = UV2.row(F2(j, 0));
            Eigen::RowVector2d B = UV2.row(F2(j, 1));
            Eigen::RowVector2d C = UV2.row(F2(j, 2));

            if (isPointInTriangle(b, A, B, C)) {
                // Register point a from V1
                pointCloudA.push_back(a);

                // Register corresponding vertices from the face in V2
                pointCloudB.push_back(V2.row(F2(j, 0)));
                pointCloudB.push_back(V2.row(F2(j, 1)));
                pointCloudB.push_back(V2.row(F2(j, 2)));

                polyscope::registerPointCloud("p1", pointCloudA);
                polyscope::registerPointCloud("p2 triangle", pointCloudB);
                

                //polyscope::show();
                //polyscope::removeAllStructures();

                break; // Move to the next vertex in V1
            }

        }
        // Progress bar logic
        int progressWidth = 50; // Width of the progress bar
        double progress = static_cast<double>(i + 1) / V1.rows();
        int pos = static_cast<int>(progress * progressWidth);

        std::cout << "\r["; // Carriage return for updating the same line
        for (int k = 0; k < progressWidth; ++k) {
            if (k < pos) std::cout << "=";
            else if (k == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << std::fixed << std::setprecision(1) << (progress * 100.0) << "% completed" << std::flush;

    }
    polyscope::show();


    // Debug output
    std::cout << "Registered point clouds:" << std::endl;
    std::cout << "PointCloudA size: " << pointCloudA.size() << std::endl;
    std::cout << "PointCloudB size: " << pointCloudB.size() << std::endl;


}
