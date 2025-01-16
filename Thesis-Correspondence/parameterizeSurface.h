#ifndef PARAMETERIZE_SURFACE_H
#define PARAMETERIZE_SURFACE_H

#include <Eigen/Dense>
#include <string>

bool parameterizeSurface(const Eigen::MatrixXd& Mesh1_V, const Eigen::MatrixXi& Mesh1_F, const std::string& output_filename);

bool parameterizeSurfaceLibiGL(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& UV);

bool parameterizeSurface3(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

bool parameterizeSurfaceTutte(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& UV);

bool parameterizeSurfaceHarmonic(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& UV);


bool paramsurface5(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& UV, const Eigen::MatrixXd& boundary_vertices, bool boundaryEnabled, const Eigen::MatrixXd& boundary_vertices_other);

void CalculateGenus(Eigen::MatrixXd V, Eigen::MatrixXi F);


Eigen::MatrixXd readAndConcatenatePointClouds(const std::string& folderPath, const std::string& filePatternStr);


void CompleteBorderCorrespondence(
    Eigen::MatrixXd& V1, Eigen::MatrixXi& F1,
    Eigen::MatrixXd& V2, Eigen::MatrixXi& F2,
    Eigen::MatrixXd& border_1, Eigen::MatrixXd& border_connected_1,
    Eigen::MatrixXd& border_2, Eigen::MatrixXd& border_connected_2
);

#endif // PARAMETERIZE_SURFACE_H
