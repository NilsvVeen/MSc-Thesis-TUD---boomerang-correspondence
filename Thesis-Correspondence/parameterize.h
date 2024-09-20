#include <Eigen/Core>
#include <Eigen/Dense>

double computeDistance(const Eigen::VectorXd& p1, const Eigen::VectorXd& p2);

Eigen::MatrixXd sortVerticesByProximity(const Eigen::MatrixXd& V);