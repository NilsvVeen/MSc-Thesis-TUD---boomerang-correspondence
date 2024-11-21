#include <iostream>
#include <vector>
#include <tuple>
#include <Eigen/Core>
#include <filesystem>
#include "stl_utils.h"
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>


void readMeshesAndPointClouds(const std::string& meshesFolder,
    const std::string& pointCloudsFolder,
    Eigen::MatrixXd& Mesh1_V, Eigen::MatrixXi& Mesh1_F,
    Eigen::MatrixXd& Mesh2_V, Eigen::MatrixXi& Mesh2_F,
    std::vector<Eigen::MatrixXd>& V1_pointclouds,
    std::vector<Eigen::MatrixXd>& V2_pointclouds);


// Function declaration
void showInPolyscope(const Eigen::MatrixXd& Mesh1_V, const Eigen::MatrixXi& Mesh1_F,
    const Eigen::MatrixXd& Mesh2_V, const Eigen::MatrixXi& Mesh2_F,
    const std::vector<Eigen::MatrixXd>& V1_pointclouds,
    const std::vector<Eigen::MatrixXd>& V2_pointclouds, bool enableMeshes = true, bool enablePointClouds = true);


void findExactCorrespondences(const Eigen::MatrixXd& mesh1_V, Eigen::MatrixXd& V1_pointcloud, std::vector<int>& indices_to_remove);


void findClosestCorrespondences(const Eigen::MatrixXd& mesh2_V, Eigen::MatrixXd& V2_pointcloud);

void writeOutputsToFolder(const std::string& outputFolder,
    const Eigen::MatrixXd& Mesh1_V, const Eigen::MatrixXi& Mesh1_F,
    const Eigen::MatrixXd& Mesh2_V, const Eigen::MatrixXi& Mesh2_F,
    const std::vector<Eigen::MatrixXd>& V1_pointclouds,
    const std::vector<Eigen::MatrixXd>& V2_pointclouds);

void projectAndSplitMesh(const Eigen::MatrixXd& mesh2_V, const Eigen::MatrixXi& mesh2_F,
    const Eigen::MatrixXd& V2_pointcloud, Eigen::MatrixXd& V2_pointcloud_new,
    Eigen::MatrixXd& updated_mesh_V, Eigen::MatrixXi& updated_mesh_F);