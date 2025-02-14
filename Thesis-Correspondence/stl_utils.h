#ifndef STL_UTILS_H
#define STL_UTILS_H

#include <iostream>
#include <filesystem>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <igl/readSTL.h>
#include <igl/writeSTL.h>
#include <igl/fit_plane.h>
#include <igl/slice_mask.h>
#include <igl/copyleft/cgal/convex_hull.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_2.h>
#include <vector>

// Function to view STL object in Polyscope
void viewSTLObject(const std::string& filename);

// Function to print files in a directory
void printFilesInDirectory(const std::string& directoryPath);

// Function to save mesh to a file in OBJ format
void saveMeshToFile(const std::string& filename, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

// Function to save point cloud to a file in PLY format
void savePointCloudToFile(const std::string& filename, const Eigen::MatrixXd& V);

// Function to initialize Polyscope and register a mesh
void initializePolyscopeAndRegisterMesh(const std::string& name, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

// Function to take a screenshot of the current Polyscope view
void takeScreenshot(const std::string& screenshotName);

// Function to show all registered meshes in Polyscope
void showMeshes();

// Function to handle STL file loading and processing
void processSTLFile(const std::string& filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& N);

Eigen::MatrixXd readVerticesFromPLY(const std::string& filename);
Eigen::MatrixXd readVerticesFromPLY2D(const std::string& filename);

void writeVerticesToPLY(const std::string& filename, const Eigen::MatrixXd& vertices);

void write2DVerticesToPLY(const std::string& filename, const Eigen::MatrixXd& vertices);

bool readMeshFromFile(const std::string& filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F);

bool readPointCloudFromFile(const std::string& filename, Eigen::MatrixXd& V);

void showPointCloudInParts(const Eigen::MatrixXd& V);

void showSelection(const Eigen::MatrixXd& V);

void showSideBySideSelection(const Eigen::MatrixXd& V1, const Eigen::MatrixXd& V2);

void showSideBySideSelectionWithVertexSelection(const Eigen::MatrixXd& V1, const Eigen::MatrixXd& V2, const std::string& dir, Eigen::Vector3d& shiftVector );

void showSideBySideMeshes(const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F1,
    const Eigen::MatrixXd& V2, const Eigen::MatrixXi& F2, const Eigen::MatrixXd V1_border, Eigen::MatrixXd V2_border, const std::string& dir);

Eigen::MatrixXd calculateAndAdjustOffsets(const Eigen::MatrixXd& V1, const Eigen::MatrixXd& V2, Eigen::MatrixXd& shiftVector);

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> calculateAndAdjustOffsetsFromBorders(
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXd& V2,
    const Eigen::MatrixXd& V1_border,
    const Eigen::MatrixXd& V2_border);

int countHoles(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

void orderGuide(Eigen::MatrixXd& V1);

void resetPolyscope();

#endif // STL_UTILS_H
