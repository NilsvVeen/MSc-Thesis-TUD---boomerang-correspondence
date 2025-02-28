// file_utils.h
#ifndef FILE_UTILS_H
#define FILE_UTILS_H

#include <vector>
#include <string>
#include "geometry.h"
#include <Eigen/Dense>

void createDirectory(const std::string& dirName);
void writeToFile(const std::string& filename, const std::vector<Point>& vertices);
void writeParamsToFile(const std::string& filename, const std::vector<double>& params);

void clearDirectory(const std::string& dirName);

void loadMeshesFromFolder(const std::string& folderPath,
    std::vector<Eigen::MatrixXd>& inputV,
    std::vector<Eigen::MatrixXi>& inputF);

void loadMeshes(const std::vector<std::string>& filePaths,
    std::vector<Eigen::MatrixXd>& inputV,
    std::vector<Eigen::MatrixXi>& inputF);


#endif // FILE_UTILS_H
