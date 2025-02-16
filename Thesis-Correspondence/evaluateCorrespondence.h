#include <Eigen/Dense>


void analyzeAndVisualizeCorrespondence(const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F1,
    const Eigen::MatrixXd& V2, const Eigen::MatrixXi& F2, const std::string evaluateCorrespondenceFolder);


void computeMeshDistances(
    const Eigen::MatrixXd& VA, const Eigen::MatrixXi& FA,
    const Eigen::MatrixXd& VB, const Eigen::MatrixXi& FB,
    double& hausdorffDist, double& chamferDist);