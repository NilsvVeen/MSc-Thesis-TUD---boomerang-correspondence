#include <Eigen/Core>
#include <Eigen/Dense>

struct Point3 {
    double x, y, z;
};

double distance2d(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2);

Eigen::MatrixXd sortVerticesByProximity(const Eigen::MatrixXd& V);

std::vector<double> computeUnitParametrization(const Eigen::MatrixXd& vertices);

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> getEqualizedPointClouds(const Eigen::MatrixXd& sortedVertices1, const Eigen::MatrixXd& sortedVertices2);