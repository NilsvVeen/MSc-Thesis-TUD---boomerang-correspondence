#ifndef PARAMETERIZE_SURFACE_H
#define PARAMETERIZE_SURFACE_H

#include <Eigen/Dense>
#include <string>

bool parameterizeSurface(const Eigen::MatrixXd& Mesh1_V, const Eigen::MatrixXi& Mesh1_F, const std::string& output_filename);

#endif // PARAMETERIZE_SURFACE_H
