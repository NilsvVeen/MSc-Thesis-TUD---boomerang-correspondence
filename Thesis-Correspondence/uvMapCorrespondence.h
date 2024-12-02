#include <Eigen/Core>       // For Eigen::MatrixXd, Eigen::MatrixXi, and basic matrix operations
#include <Eigen/Dense>      // For additional dense matrix operations (e.g., decomposition)
#include <Eigen/Geometry>   // For geometric transformations if required (e.g., rotations)



void UVToCorrespondence(
    const Eigen::MatrixXd& V1,  // Mesh 1 vertices       n x 3
    const Eigen::MatrixXi& F1, // Mesh 1 faces          m x 3
    const Eigen::MatrixXd& B1, // Boundary 1 vertices   p x 3
    const Eigen::MatrixXd& UV1, // UV map of mesh 1     n x 2

    const Eigen::MatrixXd& V2,  // Mesh 2 vertices       s x 3
    const Eigen::MatrixXi& F2,  // Mesh 2 faces          q x 3
    const Eigen::MatrixXd& B2,  // Boundary 2 vertices   p x 3
    const Eigen::MatrixXd& UV2,  // UV map of mesh 2     s x 2
    const std::string correspondence3dMatched
);
 
