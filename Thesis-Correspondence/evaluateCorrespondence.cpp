#include <Eigen/Dense>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <iostream>
#include <string>
#include <cmath>
#include "file_utils.h"


#include <fstream>
#include <numeric>




void computeAreaAndShear(const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F1,
    const Eigen::MatrixXd& V2, const Eigen::MatrixXi& F2,
    Eigen::VectorXd& areaDistortions, Eigen::VectorXd& shearDistortions) {

    int numTriangles = F1.rows();
    areaDistortions.resize(numTriangles);
    shearDistortions.resize(numTriangles);

    for (int i = 0; i < numTriangles; i++) {
        // Get the vertices of the triangle from both meshes
        Eigen::Vector3d v1_1 = V1.row(F1(i, 0));
        Eigen::Vector3d v1_2 = V1.row(F1(i, 1));
        Eigen::Vector3d v1_3 = V1.row(F1(i, 2));

        Eigen::Vector3d v2_1 = V2.row(F2(i, 0));
        Eigen::Vector3d v2_2 = V2.row(F2(i, 1));
        Eigen::Vector3d v2_3 = V2.row(F2(i, 2));

        // Compute the area of the triangle in V1 and V2
        Eigen::Vector3d edge1 = v1_2 - v1_1;
        Eigen::Vector3d edge2 = v1_3 - v1_1;
        double area1 = 0.5 * edge1.cross(edge2).norm();

        edge1 = v2_2 - v2_1;
        edge2 = v2_3 - v2_1;
        double area2 = 0.5 * edge1.cross(edge2).norm();

        // Calculate area distortion, handle degenerate triangles
        double areaDistortion = (area2 > 1e-8) ? area1 / area2 : std::numeric_limits<double>::quiet_NaN();
        areaDistortions[i] = areaDistortion;

        // Compute edge lengths in V1 and V2
        double edgeLength1_V1 = (v1_2 - v1_1).norm();
        double edgeLength2_V1 = (v1_3 - v1_2).norm();
        double edgeLength3_V1 = (v1_1 - v1_3).norm();

        double edgeLength1_V2 = (v2_2 - v2_1).norm();
        double edgeLength2_V2 = (v2_3 - v2_2).norm();
        double edgeLength3_V2 = (v2_1 - v2_3).norm();

        // Compute shear distortion
        double shearDistortion1 = (edgeLength1_V1 > 1e-8) ? edgeLength1_V2 / edgeLength1_V1 : std::numeric_limits<double>::quiet_NaN();
        double shearDistortion2 = (edgeLength2_V1 > 1e-8) ? edgeLength2_V2 / edgeLength2_V1 : std::numeric_limits<double>::quiet_NaN();
        double shearDistortion3 = (edgeLength3_V1 > 1e-8) ? edgeLength3_V2 / edgeLength3_V1 : std::numeric_limits<double>::quiet_NaN();

        // Average shear distortion, handle NaN values
        double avgShearDistortion = (std::isfinite(shearDistortion1) + std::isfinite(shearDistortion2) + std::isfinite(shearDistortion3)) > 0
            ? (std::isfinite(shearDistortion1) * shearDistortion1 +
                std::isfinite(shearDistortion2) * shearDistortion2 +
                std::isfinite(shearDistortion3) * shearDistortion3) /
            (std::isfinite(shearDistortion1) + std::isfinite(shearDistortion2) + std::isfinite(shearDistortion3))
            : std::numeric_limits<double>::quiet_NaN();

        shearDistortions[i] = avgShearDistortion;

        //std::cout << "Triangle " << i
        //    << " - Area Distortion: " << areaDistortion
        //    << ", Average Shear Distortion: " << avgShearDistortion << std::endl;
    }
}



void computeStatisticsAndHighlight(const Eigen::VectorXd& distortions,
    double nPercent,
    Eigen::VectorXd& highlightMask,
    std::string distortionType) {

    int numTriangles = distortions.size();

    // Filter valid (finite) values
    std::vector<double> validDistortions;
    for (int i = 0; i < numTriangles; ++i) {
        if (std::isfinite(distortions[i])) {
            validDistortions.push_back(distortions[i]);
        }
    }

    if (validDistortions.empty()) {
        std::cout << "No valid " << distortionType << " values to compute statistics." << std::endl;
        highlightMask.setZero(numTriangles);
        return;
    }

    // Compute statistics
    double minDistortion = *std::min_element(validDistortions.begin(), validDistortions.end());
    double maxDistortion = *std::max_element(validDistortions.begin(), validDistortions.end());
    double meanDistortion = std::accumulate(validDistortions.begin(), validDistortions.end(), 0.0) / validDistortions.size();
    double stdDevDistortion = std::sqrt(std::accumulate(validDistortions.begin(), validDistortions.end(), 0.0,
        [meanDistortion](double sum, double val) {
            return sum + (val - meanDistortion) * (val - meanDistortion);
        }) / validDistortions.size());

    // Print statistics
    std::cout << distortionType << " Distortion - Min: " << minDistortion
        << ", Max: " << maxDistortion
        << ", Mean: " << meanDistortion
        << ", StdDev: " << stdDevDistortion << std::endl;

    // Highlight top n% distortions
    std::vector<std::pair<double, int>> distortionsWithIndices;
    for (int i = 0; i < numTriangles; ++i) {
        if (std::isfinite(distortions[i])) {
            distortionsWithIndices.emplace_back(distortions[i], i);
        }
    }
    std::sort(distortionsWithIndices.rbegin(), distortionsWithIndices.rend());

    int topN = static_cast<int>(std::ceil(nPercent / 100.0 * distortionsWithIndices.size()));
    highlightMask.setZero(numTriangles);
    for (int i = 0; i < topN; ++i) {
        highlightMask[distortionsWithIndices[i].second] = 1.0; // Mark top n% triangles
    }

    std::cout << "Highlighted top " << nPercent << "% triangles with the largest "
        << distortionType << " distortion." << std::endl;
}






// Main analysis and visualization function
void analyzeAndVisualizeCorrespondence(const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F1,
    const Eigen::MatrixXd& V2, const Eigen::MatrixXi& F2, const std::string evaluateCorrespondenceFolder) {

    double nPercent = 5;

    // Compute area and shear distortions
    Eigen::VectorXd areaDistortions, shearDistortions;
    computeAreaAndShear(V1, F1, V2, F2, areaDistortions, shearDistortions);

    createDirectory(evaluateCorrespondenceFolder);
    createDirectory(evaluateCorrespondenceFolder);

    std::string areaDistortionsFile = evaluateCorrespondenceFolder + "/area_distortions.txt";
    std::ofstream areaFile(areaDistortionsFile);
    if (areaFile.is_open()) {
        areaFile << areaDistortions << "\n";
        areaFile.close();
        std::cout << "Area distortions written to: " << areaDistortionsFile << "\n";
    }
    else {
        std::cerr << "Error: Could not write to file " << areaDistortionsFile << "\n";
    }
    // Step 4: Write shear distortions to a file
    std::string shearDistortionsFile = evaluateCorrespondenceFolder + "/shear_distortions.txt";
    std::ofstream shearFile(shearDistortionsFile);
    if (shearFile.is_open()) {
        shearFile << shearDistortions << "\n";
        shearFile.close();
        std::cout << "Shear distortions written to: " << shearDistortionsFile << "\n";
    }
    else {
        std::cerr << "Error: Could not write to file " << shearDistortionsFile << "\n";
    }


    // Compute and highlight distortions separately
    Eigen::VectorXd areaHighlightMask, shearHighlightMask;
    computeStatisticsAndHighlight(areaDistortions, nPercent, areaHighlightMask, "Area");
    computeStatisticsAndHighlight(shearDistortions, nPercent, shearHighlightMask, "Shear");







    // Initialize Polyscope
    polyscope::init();

    // Register meshes
    auto* mesh1 = polyscope::registerSurfaceMesh("Mesh 1", V1, F1);
    mesh1->addFaceScalarQuantity("Area Distortion", areaDistortions);
    mesh1->addFaceScalarQuantity("Shear Distortion", shearDistortions);
    mesh1->addFaceScalarQuantity("Highlighted Area Distortion", areaHighlightMask);
    mesh1->addFaceScalarQuantity("Highlighted Shear Distortion", shearHighlightMask);
    polyscope::registerSurfaceMesh("Mesh 2", V2, F2);

    // Show the visualization
    polyscope::show();
}



#include <igl/hausdorff.h>
#include <igl/point_mesh_squared_distance.h>
#include <Eigen/Core>
#include <iostream>
#include <cmath>

// Compute the Hausdorff and Chamfer distances between two meshes:
// Mesh A: (VA, FA) and Mesh B: (VB, FB)
void computeMeshDistances(
    const Eigen::MatrixXd& VA, const Eigen::MatrixXi& FA,
    const Eigen::MatrixXd& VB, const Eigen::MatrixXi& FB,
    double& hausdorffDist, double& chamferDist)
{
    // ---------------------------------------------------------------------------
    // 1. Compute the Hausdorff distance using libigl's hausdorff function.
    //    Note: This implementation computes:
    //    d_H(A,B) = max { max_{a in A} min_{b in B} ||a-b||,  max_{b in B} min_{a in A} ||a-b|| }
    // ---------------------------------------------------------------------------
    igl::hausdorff(VA, FA, VB, FB, hausdorffDist);

    // ---------------------------------------------------------------------------
    // 2. Compute the Chamfer distance.
    //    For Chamfer, we compute the average point-to-mesh distance in both directions.
    //    We'll use igl::point_mesh_squared_distance to compute the squared distances.
    // ---------------------------------------------------------------------------
    Eigen::VectorXd sqDistA; // Squared distances from each vertex in VA to mesh (VB,FB)
    Eigen::VectorXi I_A;     // (Optional) Index of closest facet for each vertex in VA
    Eigen::MatrixXd C_A;     // (Optional) Closest point on mesh (VB,FB) for each vertex in VA

    igl::point_mesh_squared_distance(VA, VB, FB, sqDistA, I_A, C_A);

    Eigen::VectorXd sqDistB; // Squared distances from each vertex in VB to mesh (VA,FA)
    Eigen::VectorXi I_B;
    Eigen::MatrixXd C_B;

    igl::point_mesh_squared_distance(VB, VA, FA, sqDistB, I_B, C_B);

    // Convert squared distances to distances by taking square roots and sum them
    double sumA = 0.0;
    for (int i = 0; i < sqDistA.size(); ++i)
        sumA += std::sqrt(sqDistA(i));
    double sumB = 0.0;
    for (int i = 0; i < sqDistB.size(); ++i)
        sumB += std::sqrt(sqDistB(i));

    // Average the distances over the number of vertices in each mesh
    double avgA = sumA / static_cast<double>(VA.rows());
    double avgB = sumB / static_cast<double>(VB.rows());

    // Chamfer distance: average of the two directional distances.
    chamferDist = 0.5 * (avgA + avgB);



    std::cout << "Hausdorff Distance: " << hausdorffDist << std::endl;
    std::cout << "Chamfer Distance: " << chamferDist << std::endl;
}








