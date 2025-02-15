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





