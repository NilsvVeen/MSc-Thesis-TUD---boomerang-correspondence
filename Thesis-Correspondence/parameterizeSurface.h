#ifndef PARAMETERIZE_SURFACE_H
#define PARAMETERIZE_SURFACE_H

#include <Eigen/Dense>
#include <string>

bool parameterizeSurface(const Eigen::MatrixXd& Mesh1_V, const Eigen::MatrixXi& Mesh1_F, const std::string& output_filename);

bool parameterizeSurfaceLibiGL(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& UV);

bool parameterizeSurface3(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

bool parameterizeSurfaceTutte(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& UV);

bool parameterizeSurfaceHarmonic(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& UV);

bool paramsurface5(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& UV, const Eigen::MatrixXd& boundary_vertices);


Eigen::MatrixXd readAndConcatenatePointClouds(const std::string& folderPath, const std::string& filePatternStr);


#endif // PARAMETERIZE_SURFACE_H
