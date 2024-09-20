// stl_utils.h
#ifndef STL_UTILS_H
#define STL_UTILS_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <string>


void printFilesInDirectory(const std::string& directoryPath);

void viewSTLObject(const std::string& filename);

// Function to compute the surface area of a triangle
double computeTriangleArea(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2);

// Project vertices to 2D (XY plane) and fix Z to a specified value
Eigen::MatrixXd projectToXYWithFixedZ(const Eigen::MatrixXd& V, double fixedZ);

// Initialize Polyscope and register a mesh
void initializePolyscopeAndRegisterMesh(const std::string& name, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

// Take a screenshot of the current Polyscope view
void takeScreenshot(const std::string& screenshotName);

// Show all registered meshes in Polyscope
void showMeshes();

// Function to fit a plane to the mesh and align it
void fitPlaneAndAlignMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& rotatedV, Eigen::MatrixXi& rotatedF);

// Function to find the slice with maximum area
void findMaxAreaSlice(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& V_maxAreaSlice, Eigen::MatrixXi& F_maxAreaSlice, int& maxAreaSlice);

// Function to calculate average Z coordinate for a slice
double calculateAverageZ(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

// Function to handle STL file loading and processing
void processSTLFile(const std::string& filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& N);

// Main function to fit plane, align mesh, and show results
std::vector<Eigen::Vector2d> fitPlaneAndAlignMesh(const std::string& filename, const std::string& outputDir);

#endif // STL_UTILS_H
