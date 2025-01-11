#include <Polyscope/Polyscope.h>
#include <Eigen/Dense>
#include <iostream>
#include <polyscope/surface_mesh.h>

void performPCAAndEditWithVisualization(
    const std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>>& inputShapes) {

    // Ensure at least one shape is provided
    if (inputShapes.empty()) {
        std::cerr << "No input shapes provided!" << std::endl;
        return;
    }

    // Extract vertex matrices
    std::vector<Eigen::MatrixXd> vertices;
    for (const auto& shape : inputShapes) {
        vertices.push_back(shape.first);
    }

    // Compute the mean shape (assuming all shapes have the same vertex count)
    int numVertices = vertices[0].rows();
    Eigen::MatrixXd meanShape = Eigen::MatrixXd::Zero(numVertices, 3);
    for (const auto& V : vertices) {
        meanShape += V;
    }
    meanShape /= vertices.size();

    // Compute variations (centered shapes)
    Eigen::MatrixXd variations(vertices.size() * numVertices, 3);
    for (size_t i = 0; i < vertices.size(); ++i) {
        variations.block(i * numVertices, 0, numVertices, 3) = vertices[i] - meanShape;
    }

    // Flatten variations into a single matrix for PCA
    // The matrix should have shape (numVertices * numShapes) x 3
    Eigen::MatrixXd dataMatrix = variations;  // Already in correct shape (numVertices * numShapes, 3)

    // Compute PCA (Eigen decomposition of covariance matrix)
    Eigen::MatrixXd centeredData = dataMatrix.rowwise() - dataMatrix.colwise().mean();
    Eigen::MatrixXd covarianceMatrix = centeredData.transpose() * centeredData;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covarianceMatrix);

    // Eigenvalues and Eigenvectors
    Eigen::VectorXd eigenValues = eigenSolver.eigenvalues().reverse();
    Eigen::MatrixXd eigenVectors = eigenSolver.eigenvectors().rowwise().reverse();

    std::cout << "PCA computed. Eigenvalues: \n" << eigenValues.transpose() << std::endl;

    // Modes of variation
    Eigen::MatrixXd modes = eigenVectors.leftCols(3); // First 3 modes as example
    std::cout << "First 3 modes of variation computed." << std::endl;

    // Shape deformation with constraint
    Eigen::VectorXd coefficients = Eigen::VectorXd::Zero(modes.cols());
    coefficients(0) = 1.0; // Example: activate first mode
    Eigen::MatrixXd deformedShape = meanShape;
    for (int i = 0; i < modes.cols(); ++i) {
        Eigen::MatrixXd deformation = modes.col(i).replicate(numVertices, 1); // Expand mode to full vertex count
        deformedShape += coefficients(i) * deformation;
    }

    std::cout << "Deformed shape computed with constraints." << std::endl;

    // Visualize using Polyscope
    // Initialize Polyscope
    polyscope::init();

    // Add the deformed shape to Polyscope as a surface mesh
    polyscope::registerSurfaceMesh("Deformed Shape", deformedShape, inputShapes[0].second);

    // Optionally, add vertex colors (e.g., distance from mean)
    Eigen::VectorXd vertexColors(numVertices);
    for (int i = 0; i < numVertices; ++i) {
        vertexColors(i) = (deformedShape.row(i) - meanShape.row(i)).norm();
    }
    polyscope::getSurfaceMesh("Deformed Shape")->addVertexColorQuantity("Vertex Deformation", vertexColors);

    // Show the visualization
    polyscope::show();
}
