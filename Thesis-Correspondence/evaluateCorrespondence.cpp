#include <Eigen/Dense>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <iostream>
#include <string>
#include <cmath>

// Function to compute the area of a triangle
double computeTriangleArea(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int faceIndex) {
    Eigen::RowVector3d v0 = V.row(F(faceIndex, 0));
    Eigen::RowVector3d v1 = V.row(F(faceIndex, 1));
    Eigen::RowVector3d v2 = V.row(F(faceIndex, 2));

    Eigen::RowVector3d edge1 = v1 - v0;
    Eigen::RowVector3d edge2 = v2 - v0;
    double area = 0.5 * edge1.cross(edge2).norm();
    return area;
}

// Function to compute area ratios
Eigen::VectorXd computeAreaRatios(const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F1,
    const Eigen::MatrixXd& V2, const Eigen::MatrixXi& F2) {
    Eigen::VectorXd areaRatios(F1.rows());
    for (int i = 0; i < F1.rows(); ++i) {
        double area1 = computeTriangleArea(V1, F1, i);
        double area2 = computeTriangleArea(V2, F2, i);
        areaRatios(i) = area1 / area2; // Original area / Mapped area
    }
    return areaRatios;
}





void computeAreaAndShear(const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F1,
    const Eigen::MatrixXd& V2, const Eigen::MatrixXi& F2) {

    int numTriangles = F1.rows();
    double totalAreaDistortion = 0.0;
    double totalShearDistortion = 0.0;

    for (int i = 0; i < numTriangles; i++) {
        // Get the vertices of the triangle from both meshes
        Eigen::Vector3d v1_1 = V1.row(F1(i, 0));
        Eigen::Vector3d v1_2 = V1.row(F1(i, 1));
        Eigen::Vector3d v1_3 = V1.row(F1(i, 2));

        Eigen::Vector3d v2_1 = V2.row(F2(i, 0));
        Eigen::Vector3d v2_2 = V2.row(F2(i, 1));
        Eigen::Vector3d v2_3 = V2.row(F2(i, 2));

        // Compute the area of the triangle in V1
        Eigen::Vector3d edge1 = v1_2 - v1_1;
        Eigen::Vector3d edge2 = v1_3 - v1_1;
        double area1 = 0.5 * edge1.cross(edge2).norm();

        // Compute the area of the triangle in V2
        edge1 = v2_2 - v2_1;
        edge2 = v2_3 - v2_1;
        double area2 = 0.5 * edge1.cross(edge2).norm();

        // Calculate the area distortion (stretching or compression)
        double areaDistortion = area1 / area2;
        totalAreaDistortion += areaDistortion;

        // Compute edge lengths in the original mesh (V1) and transformed mesh (V2)
        double edgeLength1_V1 = edge1.norm();
        double edgeLength2_V1 = (v1_3 - v1_2).norm();
        double edgeLength3_V1 = (v1_1 - v1_3).norm();

        double edgeLength1_V2 = (v2_2 - v2_1).norm();
        double edgeLength2_V2 = (v2_3 - v2_2).norm();
        double edgeLength3_V2 = (v2_1 - v2_3).norm();

        // Compute the shear distortion based on edge lengths in the two meshes
        // Ratio of edge lengths (for each edge) could be used to compute shear distortion
        double shearDistortion1 = edgeLength1_V2 / edgeLength1_V1;
        double shearDistortion2 = edgeLength2_V2 / edgeLength2_V1;
        double shearDistortion3 = edgeLength3_V2 / edgeLength3_V1;

        // Average shear distortion for the triangle
        double avgShearDistortion = (shearDistortion1 + shearDistortion2 + shearDistortion3) / 3.0;

        totalShearDistortion += avgShearDistortion;

        std::cout << "Triangle " << i << " - Average Shear Distortion: " << avgShearDistortion << std::endl;
    }

    // Print average area and shear distortion
    double avgAreaDistortion = totalAreaDistortion / numTriangles;
    double avgShearDistortion = totalShearDistortion / numTriangles;

    std::cout << "Average Area Distortion: " << avgAreaDistortion << std::endl;
    std::cout << "Average Shear Distortion: " << avgShearDistortion << std::endl;
}




// Main analysis and visualization function
void analyzeAndVisualizeCorrespondence(const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F1,
    const Eigen::MatrixXd& V2, const Eigen::MatrixXi& F2) {

    // Compute area ratios
    Eigen::VectorXd areaRatios = computeAreaRatios(V1, F1, V2, F2);

    computeAreaAndShear(V1, F1, V2, F2);


    // Initialize Polyscope
    polyscope::init();

    // Register meshes
    polyscope::registerSurfaceMesh("Mesh 1", V1, F1)->addFaceScalarQuantity("Area Ratio", areaRatios);
    polyscope::registerSurfaceMesh("Mesh 2", V2, F2);

    // Show the visualization
    polyscope::show();
}


