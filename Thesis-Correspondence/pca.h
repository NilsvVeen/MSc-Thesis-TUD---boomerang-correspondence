#include <Eigen/Dense>
#include <iostream>

void performPCAAndEditWithVisualization(
    const std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>>& inputShapes);

std::vector<Eigen::MatrixXd> ICPAlignShapes(
    std::vector<Eigen::MatrixXd>& shapes,  // Pass shapes as non-const
    int max_iters);