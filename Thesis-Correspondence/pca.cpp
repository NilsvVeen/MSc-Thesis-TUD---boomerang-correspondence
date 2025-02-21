#include <Polyscope/Polyscope.h>
#include <Eigen/Dense>
#include <iostream>
#include <polyscope/surface_mesh.h>

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <iostream>
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/pick.h>

// Global variables for PCA results and interactive parameters:
Eigen::VectorXd g_meanShape;
Eigen::MatrixXd g_eigenvectors;
int g_numVertices = 0;
int g_vectorSize = 0;

// Interactive parameters:
int g_principalIndex = 0;
float g_weight = 0.0f;
int g_referenceShapeIndex = -1; // -1: mean shape, 0-N: input shape index

// Face connectivity and mesh storage
std::vector<std::vector<size_t>> g_faceList;
std::vector<Eigen::VectorXd> g_inputShapes; // Store input shapes as vectorized versions
polyscope::SurfaceMesh* g_deformedMesh = nullptr;

// Global variables for random shape generation
Eigen::VectorXd g_randomShape;
polyscope::SurfaceMesh* g_randomMesh = nullptr;




// Convert an Eigen vector (flattened shape) to a Polyscope vertices vector
std::vector<std::array<double, 3>> eigenVectorToVertices(const Eigen::VectorXd& shapeVec) {
    std::vector<std::array<double, 3>> vertices(g_numVertices);
    for (int v = 0; v < g_numVertices; ++v) {
        vertices[v] = { shapeVec(3 * v), shapeVec(3 * v + 1), shapeVec(3 * v + 2) };
    }
    return vertices;
}

// Generate a random shape by sampling from a Gaussian distribution
void generateRandomShape() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(-500, 500); // Standard normal

    Eigen::VectorXd randomWeights(g_eigenvectors.cols());
    for (int i = 0; i < g_eigenvectors.cols()-1; ++i) {
        randomWeights(i) = distribution(gen);
    }

    Eigen::VectorXd deformation = g_eigenvectors * randomWeights;
    g_randomShape = g_meanShape + deformation;

    // Convert to Polyscope format and update the mesh
    std::vector<std::array<double, 3>> randomVertices = eigenVectorToVertices(g_randomShape);
    if (!g_randomMesh) {
        polyscope::registerSurfaceMesh("Random Sample", randomVertices, g_faceList);
        g_randomMesh = polyscope::getSurfaceMesh("Random Sample");
    }
    else {
        g_randomMesh->updateVertexPositions(randomVertices);
    }
}

// Update the deformed mesh given the selected reference shape, principal index, and weight
void updateDeformedMesh() {
    // Select reference shape: either mean shape or one of the input shapes
    Eigen::VectorXd baseShape = (g_referenceShapeIndex == -1) ? g_meanShape : g_inputShapes[g_referenceShapeIndex];

    // Apply deformation
    Eigen::VectorXd deformation = g_eigenvectors.col(g_principalIndex) * g_weight;
    Eigen::VectorXd deformedShape = baseShape + deformation;

    // Update Polyscope mesh
    std::vector<std::array<double, 3>> deformedPos = eigenVectorToVertices(deformedShape);
    if (g_deformedMesh) {
        g_deformedMesh->updateVertexPositions(deformedPos);
    }
}

#include <polyscope/point_cloud.h>
#include "stl_utils.h"
//#include "stl_utils.h"
//
//
    void HandlUserSelectionPCA(auto& pointCloud1, const Eigen::MatrixXd & V1, std::vector<int>& selectedVertices1,
        const std::vector<std::array<double, 3>>&vertexColors1, double radius_default) {
    if (polyscope::pick::haveSelection()) {

        std::cout << "A" << std::endl;

        auto pickResult = polyscope::pick::getSelection();
        std::cout << "B" << std::endl;
        if (pickResult.first == pointCloud1 || pickResult.first->getName() == "Point Cloud") {
            std::cout << "C" << std::endl;

            handleSelection(pointCloud1, selectedVertices1, vertexColors1);
        }

    }

    std::cout << "Selected vertices in Point Cloud 1: ";
    for (int v : selectedVertices1) std::cout << v << " ";
    std::cout << std::endl;
}




void generateAndVisualizePoints(const Eigen::MatrixXd& M) {
    // Get min and max X, Y coordinates
    double minX = M.col(0).minCoeff();
    double maxX = M.col(0).maxCoeff();
    double minY = M.col(1).minCoeff();
    double maxY = M.col(1).maxCoeff();

    // Calculate average Z value
    double avgZ = M.col(2).mean();

    // Print the min/max values and the average Z for debugging
    std::cout << "Min X: " << minX << ", Max X: " << maxX << "\n";
    std::cout << "Min Y: " << minY << ", Max Y: " << maxY << "\n";
    std::cout << "Average Z: " << avgZ << "\n";

    // Make the bounding box bigger by expanding it
    double dx = (maxX - minX) * 1.0;  
    double dy = (maxY - minY) * 1.0;  

    // Expand the bounding box in all directions
    minX -= dx;
    maxX += dx;
    minY -= dy;
    maxY += dy;

    // Print the updated min/max values for debugging
    std::cout << "Expanded Min X: " << minX << ", Expanded Max X: " << maxX << "\n";
    std::cout << "Expanded Min Y: " << minY << ", Expanded Max Y: " << maxY << "\n";

    // Number of points to generate (total number of points)
    int numPoints = 5000;

    // We will divide the expanded bounding box into a grid of approximately sqrt(numPoints) points per dimension
    int gridSize = std::sqrt(numPoints);

    // Step sizes for X and Y dimensions
    double stepX = (maxX - minX) / (gridSize - 1);
    double stepY = (maxY - minY) / (gridSize - 1);

    // Generate points in the expanded X-Y bounding box with the average Z value
    std::vector<Eigen::Vector3d> points;
    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            // Calculate X and Y coordinates based on the grid
            double x = minX + i * stepX;
            double y = minY + j * stepY;

            // Use average Z value for all points
            Eigen::Vector3d point(x, y, avgZ);

            // Print the generated point
            std::cout << "Generated Point: (" << point[0] << ", " << point[1] << ", " << point[2] << ")\n";

            // Store the point
            points.push_back(point);
        }
    }

    // Prepare points as Eigen::MatrixXd
    Eigen::MatrixXd polyscopePoints(points.size(), 3);
    for (size_t i = 0; i < points.size(); ++i) {
        polyscopePoints(i, 0) = points[i][0];
        polyscopePoints(i, 1) = points[i][1];
        polyscopePoints(i, 2) = points[i][2];
    }

    std::vector<int> selectedVerticesXXX = { 1 };

    std::cout << "selectv:" << selectedVerticesXXX[0] << " " << selectedVerticesXXX.size() << std::endl;

    std::vector<std::array<double, 3>> vertexColors1  = std::vector<std::array<double, 3>>(polyscopePoints.rows(), { {1.0, 1.0, 1.0} });
    ;
    double radius = 0.005;
    // Visualize the points in Polyscope
    //auto obj1 = polyscope::registerPointCloud("Generated Points", polyscopePoints);
    auto* obj1 = registerPointCloudWithColors("Point Cloud", polyscopePoints, radius, vertexColors1);
    obj1->updatePointPositions(polyscopePoints); // Update Polyscope 
    polyscope::state::userCallback = [&obj1, &polyscopePoints, &selectedVerticesXXX, &vertexColors1, &radius]() {

        std::cout << "selectv:" << selectedVerticesXXX[0] << " " << selectedVerticesXXX.size() << std::endl;

        HandlUserSelectionPCA(obj1, polyscopePoints, selectedVerticesXXX, vertexColors1, radius);
 
    };
    polyscope::show();



}


// Main PCA computation and visualization setup
void performPCAAndEditWithVisualization(const std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>>& inputShapes) {
    if (inputShapes.empty()) {
        std::cerr << "No shapes provided." << std::endl;
        return;
    }

    int numShapes = static_cast<int>(inputShapes.size());
    g_numVertices = static_cast<int>(inputShapes[0].first.rows());
    g_vectorSize = 3 * g_numVertices;

    // Store input shapes
    Eigen::MatrixXd data(g_vectorSize, numShapes);
    g_inputShapes.clear();
    for (int i = 0; i < numShapes; ++i) {
        Eigen::VectorXd shapeVec(g_vectorSize);
        for (int v = 0; v < g_numVertices; ++v) {
            shapeVec.segment<3>(3 * v) = inputShapes[i].first.row(v).transpose();
        }
        data.col(i) = shapeVec;
        g_inputShapes.push_back(shapeVec);
    }

    // Compute mean shape and center data
    g_meanShape = data.rowwise().mean();
    Eigen::MatrixXd centered = data.colwise() - g_meanShape;

    // Compute PCA using SVD
    Eigen::BDCSVD<Eigen::MatrixXd> svd(centered, Eigen::ComputeThinU);
    g_eigenvectors = svd.matrixU();

    std::cout << "SVD complete. Principal modes: " << g_eigenvectors.cols() << std::endl;

    // Default selection
    g_principalIndex = g_eigenvectors.cols() - 1;
    g_weight = 0.0f;

    // Convert faces to Polyscope format
    g_faceList.clear();
    for (int i = 0; i < inputShapes[0].second.rows(); ++i) {
        std::vector<size_t> face(inputShapes[0].second.cols());
        for (int j = 0; j < inputShapes[0].second.cols(); ++j) {
            face[j] = static_cast<size_t>(inputShapes[0].second(i, j));
        }
        g_faceList.push_back(face);
    }

    // Initialize Polyscope
    polyscope::init();

    // Register the mean shape
    polyscope::registerSurfaceMesh("Mean Shape", eigenVectorToVertices(g_meanShape), g_faceList);

    // Register first deformed shape
    polyscope::registerSurfaceMesh("Deformed Shape", eigenVectorToVertices(g_meanShape), g_faceList);
    g_deformedMesh = polyscope::getSurfaceMesh("Deformed Shape");

    // User controls via ImGui
    polyscope::state::userCallback = []() {
        int totalModes = g_eigenvectors.cols();
        // Dropdown to choose reference shape
        if (ImGui::BeginCombo("Reference Shape", g_referenceShapeIndex == -1 ? "Mean Shape" : ("Input Shape " + std::to_string(g_referenceShapeIndex)).c_str())) {
            if (ImGui::Selectable("Mean Shape", g_referenceShapeIndex == -1)) {
                g_referenceShapeIndex = -1;
            }
            for (int i = 0; i < g_inputShapes.size(); ++i) {
                if (ImGui::Selectable(("Input Shape " + std::to_string(i)).c_str(), g_referenceShapeIndex == i)) {
                    g_referenceShapeIndex = i;
                }
            }
            ImGui::EndCombo();
        }
        // Slider for mode selection
        ImGui::SliderInt("Principal Mode", &g_principalIndex, 0, totalModes - 1);
        // Slider for weight
        ImGui::SliderFloat("Weight", &g_weight, -10000.0f, 10000.0f);
        // Reset weight button
        if (ImGui::Button("Reset Weight")) {
            g_weight = 0.0f;
        }
        // Generate random sample button
        if (ImGui::Button("Generate Random Sample")) {
            generateRandomShape();
        }
        // Update deformed shape
        updateDeformedMesh();
    };



    //generateAndVisualizePoints(g_meanShape);
    //generateAndVisualizePoints(inputShapes[0].first);

    // Show Polyscope UI
    polyscope::show();
}






//void performPCAAndEditWithVisualization(
//    const std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>>& inputShapes) {
//
//    // Ensure at least one shape is provided
//    if (inputShapes.empty()) {
//        std::cerr << "No input shapes provided!" << std::endl;
//        return;
//    }
//
//    // Extract vertex matrices
//    std::vector<Eigen::MatrixXd> vertices;
//    for (const auto& shape : inputShapes) {
//        vertices.push_back(shape.first);
//    }
//
//    // Compute the mean shape (assuming all shapes have the same vertex count)
//    int numVertices = vertices[0].rows();
//    Eigen::MatrixXd meanShape = Eigen::MatrixXd::Zero(numVertices, 3);
//    for (const auto& V : vertices) {
//        meanShape += V;
//    }
//    meanShape /= vertices.size();
//
//    // Compute variations (centered shapes)
//    Eigen::MatrixXd variations(vertices.size() * numVertices, 3);
//    for (size_t i = 0; i < vertices.size(); ++i) {
//        variations.block(i * numVertices, 0, numVertices, 3) = vertices[i] - meanShape;
//    }
//
//    // Flatten variations into a single matrix for PCA
//    // The matrix should have shape (numVertices * numShapes) x 3
//    Eigen::MatrixXd dataMatrix = variations;  // Already in correct shape (numVertices * numShapes, 3)
//
//    // Compute PCA (Eigen decomposition of covariance matrix)
//    Eigen::MatrixXd centeredData = dataMatrix.rowwise() - dataMatrix.colwise().mean();
//    Eigen::MatrixXd covarianceMatrix = centeredData.transpose() * centeredData;
//    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covarianceMatrix);
//
//    // Eigenvalues and Eigenvectors
//    Eigen::VectorXd eigenValues = eigenSolver.eigenvalues().reverse();
//    Eigen::MatrixXd eigenVectors = eigenSolver.eigenvectors().rowwise().reverse();
//
//    std::cout << "PCA computed. Eigenvalues: \n" << eigenValues.transpose() << std::endl;
//
//    // Modes of variation
//    Eigen::MatrixXd modes = eigenVectors.leftCols(3); // First 3 modes as example
//    std::cout << "First 3 modes of variation computed." << std::endl;
//
//    // Shape deformation with constraint
//    Eigen::VectorXd coefficients = Eigen::VectorXd::Zero(modes.cols());
//    coefficients(0) = 1.0; // Example: activate first mode
//    Eigen::MatrixXd deformedShape = meanShape;
//    for (int i = 0; i < modes.cols(); ++i) {
//        Eigen::MatrixXd deformation = modes.col(i).replicate(numVertices, 1); // Expand mode to full vertex count
//        deformedShape += coefficients(i) * deformation;
//    }
//
//    std::cout << "Deformed shape computed with constraints." << std::endl;
//
//    // Visualize using Polyscope
//    // Initialize Polyscope
//    polyscope::init();
//
//    // Add the deformed shape to Polyscope as a surface mesh
//    polyscope::registerSurfaceMesh("Deformed Shape", deformedShape, inputShapes[0].second);
//
//    // Optionally, add vertex colors (e.g., distance from mean)
//    Eigen::VectorXd vertexColors(numVertices);
//    for (int i = 0; i < numVertices; ++i) {
//        vertexColors(i) = (deformedShape.row(i) - meanShape.row(i)).norm();
//    }
//    polyscope::getSurfaceMesh("Deformed Shape")->addVertexColorQuantity("Vertex Deformation", vertexColors);
//
//    // Show the visualization
//    polyscope::show();
//}



#include <Eigen/Dense>
#include <vector>
#include <limits>

// Function to find closest points in target for each point in source
std::vector<int> findClosestPoints(const Eigen::MatrixXd& source, const Eigen::MatrixXd& target) {
    std::vector<int> closestIndices(source.rows());

    for (int i = 0; i < source.rows(); ++i) {
        double minDist = std::numeric_limits<double>::max();
        int bestIndex = -1;

        for (int j = 0; j < target.rows(); ++j) {
            double dist = (source.row(i) - target.row(j)).squaredNorm();
            if (dist < minDist) {
                minDist = dist;
                bestIndex = j;
            }
        }

        closestIndices[i] = bestIndex;
    }
    return closestIndices;
}

// Function to compute optimal rigid transformation
void computeRigidTransform(const Eigen::MatrixXd& source, const Eigen::MatrixXd& target, Eigen::Matrix3d& R, Eigen::Vector3d& t) {
    Eigen::Vector3d centroidSource = source.colwise().mean();
    Eigen::Vector3d centroidTarget = target.colwise().mean();

    Eigen::MatrixXd centeredSource = source.rowwise() - centroidSource.transpose();
    Eigen::MatrixXd centeredTarget = target.rowwise() - centroidTarget.transpose();

    Eigen::Matrix3d H = centeredSource.transpose() * centeredTarget;
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);

    R = svd.matrixV() * svd.matrixU().transpose();
    if (R.determinant() < 0) {
        Eigen::Matrix3d V = svd.matrixV();
        V.col(2) *= -1;
        R = V * svd.matrixU().transpose();
    }

    t = centroidTarget - R * centroidSource;
}

// ICP function
std::vector<Eigen::MatrixXd> ICPAlignShapes(
    std::vector<Eigen::MatrixXd>& shapes,  // Pass shapes as non-const
    int max_iters) {

    if (shapes.empty()) return {};

    Eigen::MatrixXd reference = shapes[0]; // Use first shape as reference
    std::vector<Eigen::MatrixXd> aligned_shapes(shapes.size());
    aligned_shapes[0] = reference; // First shape remains unchanged

    for (size_t i = 1; i < shapes.size(); ++i) {
        Eigen::MatrixXd source = shapes[i];

        for (int iter = 0; iter < max_iters; ++iter) {
            std::vector<int> closestIndices = findClosestPoints(source, reference);

            Eigen::MatrixXd correspondingPoints(source.rows(), 3);
            for (size_t j = 0; j < closestIndices.size(); ++j) {
                correspondingPoints.row(j) = reference.row(closestIndices[j]);
            }

            Eigen::Matrix3d R;
            Eigen::Vector3d t;
            computeRigidTransform(source, correspondingPoints, R, t);

            source = (source * R.transpose()).rowwise() + t.transpose();
        }

        aligned_shapes[i] = source;
    }
    return aligned_shapes;
}





