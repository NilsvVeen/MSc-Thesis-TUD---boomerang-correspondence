#include <Polyscope/Polyscope.h>
#include <Eigen/Dense>
#include <iostream>
#include <polyscope/surface_mesh.h>

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// Function to perform PCA and prepare for visualization/editing.
void performPCAAndEditWithVisualization(const std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>>& inputShapes) {
    if (inputShapes.empty()) {
        std::cerr << "No shapes provided." << std::endl;
        return;
    }

    // Assuming all shapes have the same number of vertices.
    const int numShapes = static_cast<int>(inputShapes.size());
    const int numVertices = static_cast<int>(inputShapes[0].first.rows());
    const int vectorSize = 3 * numVertices; // each shape becomes a vector of length 3*numVertices

    // Build data matrix where each column is the vectorized shape.
    Eigen::MatrixXd data(vectorSize, numShapes);

    for (int i = 0; i < numShapes; ++i) {
        const Eigen::MatrixXd& vertices = inputShapes[i].first; // vertices: numVertices x 3
        Eigen::VectorXd shapeVec(vectorSize);
        for (int v = 0; v < numVertices; ++v) {
            // Place vertex v (x, y, z) into the vector.
            shapeVec.segment<3>(3 * v) = vertices.row(v).transpose();
        }
        data.col(i) = shapeVec;
    }



    // Compute the mean shape.
    Eigen::VectorXd meanShape = data.rowwise().mean();

    // Center the data by subtracting the mean shape from each shape.
    Eigen::MatrixXd centered = data.colwise() - meanShape;

    // Compute the covariance matrix.
    // Note: using (n-1) for an unbiased estimate.



    // Compute SVD on the centered data. We compute thin SVD for efficiency.
    // Note: The SVD of X (centered data) is X = U * S * V^T.
    // The columns of U are the principal directions (eigenvectors of the covariance matrix),
    // and the squared singular values (S.array().square()) are proportional to the variance
    // along each principal direction.
// Compute SVD on the centered data (each column is a centered shape)
    Eigen::BDCSVD<Eigen::MatrixXd> svd(centered, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Retrieve the singular values and principal components (eigenvectors of the covariance matrix)
    Eigen::VectorXd singularValues = svd.singularValues();
    Eigen::MatrixXd eigenvectors = svd.matrixU();

    // Compute the eigenvalues of the covariance matrix
    // They are given by (singularValue)^2 / (numShapes - 1)
    Eigen::VectorXd eigenvalues = singularValues.array().square() / (numShapes - 1);

    // Output the eigenvalues and eigenvectors
    //std::cout << "Eigenvalues:\n" << eigenvalues << "\n\n";
    //std::cout << "Eigenvectors (columns):\n" << eigenvectors << std::endl;



    std::cout << "SVD complete" << std::endl;



    // Use the number of columns in the eigenvectors matrix for the principal index.
    int principalIndex = eigenvectors.cols() - 1;  // Correct indexing

    // Example: create a deformed shape by moving one standard deviation along the principal component.
    double weight = 1.0;  // adjust the weight as needed
    Eigen::VectorXd deformation = eigenvectors.col(principalIndex) * weight;
    Eigen::VectorXd deformedShape = meanShape + deformation;

    // Convert the deformed shape vector back to a vertices matrix.
    Eigen::MatrixXd deformedVertices(numVertices, 3);
    for (int v = 0; v < numVertices; ++v) {
        deformedVertices.row(v) = deformedShape.segment<3>(3 * v).transpose();
    }

    // Visualization (pseudocode):
    // visualizeShapes(meanShape, deformedVertices, inputShapes[0].second);
    // Here, you would use your preferred visualization library to display the mean shape,
    // the deformed shape, and perhaps allow user interaction to adjust the weight along principal modes.

    std::cout << "PCA and deformation computed. Ready for visualization/editing." << std::endl;

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





