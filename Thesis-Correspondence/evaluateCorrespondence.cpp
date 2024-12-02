#include <Eigen/Dense>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <iostream>
#include <string>
#include <cmath>







// Function to compute area and shear distortions
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

        // Calculate area distortion
        double areaDistortion = area1 / area2;
        areaDistortions[i] = areaDistortion;

        // Compute edge lengths in V1 and V2
        double edgeLength1_V1 = (v1_2 - v1_1).norm();
        double edgeLength2_V1 = (v1_3 - v1_2).norm();
        double edgeLength3_V1 = (v1_1 - v1_3).norm();

        double edgeLength1_V2 = (v2_2 - v2_1).norm();
        double edgeLength2_V2 = (v2_3 - v2_2).norm();
        double edgeLength3_V2 = (v2_1 - v2_3).norm();

        // Compute shear distortion
        double shearDistortion1 = edgeLength1_V2 / edgeLength1_V1;
        double shearDistortion2 = edgeLength2_V2 / edgeLength2_V1;
        double shearDistortion3 = edgeLength3_V2 / edgeLength3_V1;

        // Average shear distortion
        double avgShearDistortion = (shearDistortion1 + shearDistortion2 + shearDistortion3) / 3.0;
        shearDistortions[i] = avgShearDistortion;

        std::cout << "Triangle " << i
            << " - Area Distortion: " << areaDistortion
            << ", Average Shear Distortion: " << avgShearDistortion << std::endl;
    }

    // Compute and print average metrics
    double avgAreaDistortion = areaDistortions.mean();
    double avgShearDistortion = shearDistortions.mean();

    std::cout << "Average Area Distortion: " << avgAreaDistortion << std::endl;
    std::cout << "Average Shear Distortion: " << avgShearDistortion << std::endl;
}


// Reusable function to compute statistics and highlight top n% distortions
void computeStatisticsAndHighlight(const Eigen::VectorXd& distortions,
    double nPercent,
    Eigen::VectorXd& highlightMask,
    std::string distortionType) {
    int numTriangles = distortions.size();
    int topN = static_cast<int>(std::ceil(nPercent / 100.0 * numTriangles));

    // Compute statistics
    double minDistortion = distortions.minCoeff();
    double maxDistortion = distortions.maxCoeff();
    double meanDistortion = distortions.mean();
    double stdDevDistortion = std::sqrt((distortions.array() - meanDistortion).square().mean());

    // Print statistics
    std::cout << distortionType << " Distortion - Min: " << minDistortion
        << ", Max: " << maxDistortion
        << ", Mean: " << meanDistortion
        << ", StdDev: " << stdDevDistortion << std::endl;

    // Highlight top n% distortions
    std::vector<std::pair<double, int>> distortionsWithIndices;
    for (int i = 0; i < numTriangles; ++i) {
        distortionsWithIndices.emplace_back(distortions[i], i);
    }
    std::sort(distortionsWithIndices.rbegin(), distortionsWithIndices.rend());

    highlightMask.setZero(numTriangles);
    for (int i = 0; i < topN; ++i) {
        highlightMask[distortionsWithIndices[i].second] = 1.0; // Mark top n% triangles
    }

    std::cout << "Highlighted top " << nPercent << "% triangles with the largest "
        << distortionType << " distortion." << std::endl;
}




// Main analysis and visualization function
void analyzeAndVisualizeCorrespondence(const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F1,
    const Eigen::MatrixXd& V2, const Eigen::MatrixXi& F2) {

    double nPercent = 5;

    // Compute area and shear distortions
    Eigen::VectorXd areaDistortions, shearDistortions;
    computeAreaAndShear(V1, F1, V2, F2, areaDistortions, shearDistortions);

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


