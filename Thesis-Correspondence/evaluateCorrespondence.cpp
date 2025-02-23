#include <Eigen/Dense>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <iostream>
#include <string>
#include <cmath>
#include "file_utils.h"


#include <fstream>
#include <numeric>

double computeSurfaceArea(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    double totalArea = 0.0;

    for (int i = 0; i < F.rows(); ++i) {
        Eigen::Vector3d v0 = V.row(F(i, 0));
        Eigen::Vector3d v1 = V.row(F(i, 1));
        Eigen::Vector3d v2 = V.row(F(i, 2));

        // Compute the cross product of the two edge vectors
        Eigen::Vector3d crossProduct = (v1 - v0).cross(v2 - v0);

        // Triangle area = 0.5 * norm of cross product
        totalArea += 0.5 * crossProduct.norm();
    }

    return totalArea;
}


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



std::vector<double> computeStatisticsAndHighlight(const Eigen::VectorXd& distortions,
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
        return {};
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


    return std::vector<double>{minDistortion, maxDistortion, meanDistortion, stdDevDistortion};

}






// Main analysis and visualization function
void analyzeAndVisualizeCorrespondence(const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F1,
    const Eigen::MatrixXd& V2, const Eigen::MatrixXi& F2, const std::string evaluateCorrespondenceFolder) {

    double nPercent = 5;

    // Compute area and shear distortions
    Eigen::VectorXd areaDistortions, shearDistortions;
    computeAreaAndShear(V1, F1, V2, F2, areaDistortions, shearDistortions);

    double area1 = computeSurfaceArea(V1,F1);
    double area2 = computeSurfaceArea(V2,F2);



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
    std::vector<double> areaStats = computeStatisticsAndHighlight(areaDistortions, nPercent, areaHighlightMask, "Area");
    std::vector<double> shearStats = computeStatisticsAndHighlight(shearDistortions, nPercent, shearHighlightMask, "Shear");


    // Step 4: Write shear distortions to a file
    std::string area_shear_summaryFile = evaluateCorrespondenceFolder + "/area_shear_summary.txt";
    std::ofstream sumFile(area_shear_summaryFile);
    if (sumFile.is_open()) {

        sumFile << "Area Distortion" << "\n";
        sumFile << "Distortion - Min: " << areaStats[0] << "\n";
        sumFile << "Max: " << areaStats[1] << "\n";
        sumFile << "Mean: " << areaStats[2] << "\n";
        sumFile << "StdDev " << areaStats[3] << "\n";
        

        sumFile << "Shear Distortion" << "\n";
        sumFile << "Distortion - Min: " << shearStats[0] << "\n";
        sumFile << "Max: " << shearStats[1] << "\n";
        sumFile << "Mean: " << shearStats[2] << "\n";
        sumFile << "StdDev " << shearStats[3] << "\n";

        sumFile << "Area Object 1:" << area1 << "\n";
        sumFile << "Area Object 2:" << area2 << "\n";
        sumFile.close();


        std::cout << "Distortions written to: " << area_shear_summaryFile << "\n";
    }
    else {
        std::cerr << "Error: Could not write to file " << area_shear_summaryFile << "\n";
    }




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
#include <igl/bounding_box_diagonal.h>

// Compute the Hausdorff and Chamfer distances between two meshes:
// Mesh A: (VA, FA) and Mesh B: (VB, FB)
void computeMeshDistances(
    const Eigen::MatrixXd& VA, const Eigen::MatrixXi& FA,
    const Eigen::MatrixXd& VB, const Eigen::MatrixXi& FB,
    double& hausdorffDist, double& chamferDist, const std::string evaluateCorrespondenceFolder)
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



    // Convert squared distances to actual distances.
    Eigen::VectorXd errors = sqDistB.array().sqrt();

    // ---------------------------------------------------------------------------
    // Compute extra statistics:
    //   - Mean, median, minimum, maximum, and standard deviation of the errors.
    double meanError = errors.mean();

    // For median, copy errors into a std::vector and sort.
    std::vector<double> errorVec(errors.data(), errors.data() + errors.size());
    std::sort(errorVec.begin(), errorVec.end());
    double medianError = errorVec[errorVec.size() / 2];
    double minError = errorVec.front();
    double maxError = errorVec.back();

    double sumSq = 0.0;
    for (int i = 0; i < errors.size(); i++)
        sumSq += std::pow(errors(i) - meanError, 2);
    double stdDevError = std::sqrt(sumSq / errors.size());

    // Compute bounding box diagonal of the original mesh (VA) to normalize errors.
    double bboxDiag = igl::bounding_box_diagonal(VA);

    std::cout << "Error Statistics (in original mesh units):\n";
    std::cout << "  Mean Error: " << meanError << "\n";
    std::cout << "  Median Error: " << medianError << "\n";
    std::cout << "  Min Error: " << minError << "\n";
    std::cout << "  Max Error: " << maxError << "\n";
    std::cout << "  Standard Deviation: " << stdDevError << "\n\n";

    std::cout << "Normalized by bounding box diagonal (" << bboxDiag << "):\n";
    std::cout << "  Normalized Hausdorff Distance: " << hausdorffDist / bboxDiag << "\n";
    std::cout << "  Normalized Chamfer Distance: " << chamferDist / bboxDiag << "\n";
    std::cout << "  Normalized Mean Error: " << meanError / bboxDiag << "\n";
    std::cout << "  Normalized Median Error: " << medianError / bboxDiag << "\n";
    std::cout << "  Normalized Min Error: " << minError / bboxDiag << "\n";
    std::cout << "  Normalized Max Error: " << maxError / bboxDiag << "\n";
    std::cout << "  Normalized Std Dev: " << stdDevError / bboxDiag << "\n";


    // Step 4: Write shear distortions to a file
    std::string hauschamSummaryFile = evaluateCorrespondenceFolder + "/HausChamDistancesSummary.txt";
    std::ofstream sumFile(hauschamSummaryFile);
    if (sumFile.is_open()) {

        sumFile << "Hausdorff Distance: " << hausdorffDist << "\n";
        sumFile << "Chamfer Distance: " << chamferDist << "\n";
        sumFile << "Normalized Hausdorff Distance: " << hausdorffDist / bboxDiag << "\n";
        sumFile << "Normalized Chamfer Distance: " << chamferDist / bboxDiag << "\n";

        sumFile.close();


        std::cout << "Distances written to: " << hauschamSummaryFile << "\n";
    }
    else {
        std::cerr << "Error: Could not write to file " << hauschamSummaryFile << "\n";
    }



    // ---------------------------------------------------------------------------
    // Visualize the error distribution on the deformed mesh (VB) using Polyscope.
    // We color the mesh vertices based on their error value.
    polyscope::init();


    // Register the deformed mesh in Polyscope.
    auto* psMesh = polyscope::registerSurfaceMesh("Deformed Mesh", VB, FB);

    // Optionally, you can also add a scalar quantity (the error values) so you see a legend.
    psMesh->addVertexScalarQuantity("Error Values", errors);

    polyscope::show();








}











