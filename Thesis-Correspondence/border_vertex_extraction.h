#ifndef MESH_PROCESSING_H
#define MESH_PROCESSING_H

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



// Function to compute the surface area of a triangle
double computeTriangleArea(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2);

// Project vertices to 2D (XY plane) and fix Z to a specified value
Eigen::MatrixXd projectToXYWithFixedZ(const Eigen::MatrixXd& V, double fixedZ);

// Function to fit a plane to the mesh and align it
void fitPlaneAndAlignMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& rotatedV, Eigen::MatrixXi& rotatedF);

// Function to find the slice with maximum area
void findMaxAreaSlice(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& V_maxAreaSlice, Eigen::MatrixXi& F_maxAreaSlice, int& maxAreaSlice);

// Function to calculate average Z coordinate for a slice
double calculateAverageZ(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

// Function to find border vertices using Alpha Shape
void findBorderVerticesWithAlphaShape(const Eigen::MatrixXd& V_2D, std::vector<Eigen::Vector2d>& borderVertices, float& alpha);

// Main function to fit plane, align mesh, and show results
//std::vector<Eigen::Vector2d> fitPlaneAndAlignMesh(const std::string& filename, const std::string& outputDir);
std::vector<Eigen::Vector2d> fitPlaneAndAlignMesh(const std::string& filename, const std::string& outputDir, float& alpha, Eigen::MatrixXd V_other = Eigen::MatrixXd::Zero(3, 3), Eigen::MatrixXd B_other = Eigen::MatrixXd::Zero(3, 3), bool shift = false);


#endif // MESH_PROCESSING_H
