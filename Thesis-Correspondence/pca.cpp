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

// Global variables for the PCA result and interactive parameters:
Eigen::VectorXd g_meanShape;
Eigen::MatrixXd g_eigenvectors;
int g_numVertices = 0;
int g_vectorSize = 0;

// Interactive parameters for mode selection and weight:
int g_principalIndex = 0;
float g_weight = 0.0f;

// The face connectivity of the mesh (assumed same for all shapes).
std::vector<std::vector<size_t>> g_faceList;

// Pointer to the Polyscope surface mesh so we can update its vertices.
polyscope::SurfaceMesh* g_deformedMesh = nullptr;

// Function to convert an Eigen vector (flattened shape) to a vertices vector for Polyscope.
std::vector<std::array<double, 3>> eigenVectorToVertices(const Eigen::VectorXd& shapeVec) {
    std::vector<std::array<double, 3>> vertices;
    vertices.reserve(g_numVertices);
    for (int v = 0; v < g_numVertices; ++v) {
        vertices.push_back({ shapeVec(3 * v), shapeVec(3 * v + 1), shapeVec(3 * v + 2) });
    }
    return vertices;
}

// Function to update the deformed mesh given the current g_weight and g_principalIndex.
void updateDeformedMesh() {
    // Create the deformed shape by moving along the chosen principal mode.
    Eigen::VectorXd deformation = g_eigenvectors.col(g_principalIndex) * g_weight;
    Eigen::VectorXd deformedShape = g_meanShape + deformation;

    // Convert deformedShape to vertices and update the mesh.
    std::vector<std::array<double, 3>> deformedPos = eigenVectorToVertices(deformedShape);
    if (g_deformedMesh) {
        g_deformedMesh->updateVertexPositions(deformedPos);
    }
}

// The main function that computes PCA and sets up visualization.
void performPCAAndEditWithVisualization(const std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>>& inputShapes) {
    if (inputShapes.empty()) {
        std::cerr << "No shapes provided." << std::endl;
        return;
    }

    const int numShapes = static_cast<int>(inputShapes.size());
    g_numVertices = static_cast<int>(inputShapes[0].first.rows());
    g_vectorSize = 3 * g_numVertices; // each shape becomes a vector of length 3*numVertices

    // Build data matrix where each column is the vectorized shape.
    Eigen::MatrixXd data(g_vectorSize, numShapes);
    for (int i = 0; i < numShapes; ++i) {
        const Eigen::MatrixXd& vertices = inputShapes[i].first; // vertices: numVertices x 3
        Eigen::VectorXd shapeVec(g_vectorSize);
        for (int v = 0; v < g_numVertices; ++v) {
            // Place vertex v (x, y, z) into the vector.
            shapeVec.segment<3>(3 * v) = vertices.row(v).transpose();
        }
        data.col(i) = shapeVec;
    }

    // Compute the mean shape.
    g_meanShape = data.rowwise().mean();

    // Center the data.
    Eigen::MatrixXd centered = data.colwise() - g_meanShape;

    // Instead of computing the full covariance matrix (which is huge), we compute an SVD on the centered data.
    // Using BDCSVD for efficiency (especially for tall matrices) and computing the thin SVD.
    Eigen::BDCSVD<Eigen::MatrixXd> svd(centered, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // The columns of U (stored in matrixU()) are the principal directions.
    g_eigenvectors = svd.matrixU();  // Each column is a principal direction.
    // You can also compute the eigenvalues if desired:
    // Eigen::VectorXd singularValues = svd.singularValues();
    // Eigen::VectorXd eigenvalues = singularValues.array().square() / (numShapes - 1);

    std::cout << "SVD complete. Number of principal modes: " << g_eigenvectors.cols() << std::endl;

    // Set default principal index to the last (largest variance) mode.
    g_principalIndex = static_cast<int>(g_eigenvectors.cols()) - 1;
    g_weight = 0.0f;

    // Convert the mean shape to a vertices vector.
    std::vector<std::array<double, 3>> meanPos = eigenVectorToVertices(g_meanShape);

    // Convert face indices from Eigen::MatrixXi to vector of vectors.
    const Eigen::MatrixXi& facesMat = inputShapes[0].second;
    g_faceList.clear();
    for (int i = 0; i < facesMat.rows(); ++i) {
        std::vector<size_t> face;
        for (int j = 0; j < facesMat.cols(); ++j) {
            face.push_back(static_cast<size_t>(facesMat(i, j)));
        }
        g_faceList.push_back(face);
    }

    // Initialize Polyscope.
    polyscope::init();

    // Register the mean shape for reference.
    polyscope::registerSurfaceMesh("Mean Shape", meanPos, g_faceList);

    // Compute initial deformed shape (with zero weight) and register.
    std::vector<std::array<double, 3>> deformedPos = eigenVectorToVertices(g_meanShape);
    polyscope::registerSurfaceMesh("Deformed Shape", deformedPos, g_faceList);
    g_deformedMesh = polyscope::getSurfaceMesh("Deformed Shape");

    // Add a user callback (ImGui) to control PCA modes.
    polyscope::state::userCallback = []() {
        // Get the total number of modes.
        int totalModes = g_eigenvectors.cols();

        // Slider to choose the principal mode index.
        ImGui::SliderInt("Principal Mode", &g_principalIndex, 0, totalModes - 1);

        // Slider to adjust the weight (amount of variation).
        ImGui::SliderFloat("Weight", &g_weight, -10000.0f, 10000.0f);

        // Button to reset weight.
        if (ImGui::Button("Reset Weight")) {
            g_weight = 0.0f;
        }

        // Update the deformed mesh based on current settings.
        updateDeformedMesh();
    };

    // Show the Polyscope GUI.
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





