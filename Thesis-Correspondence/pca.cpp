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
std::vector<float> g_weights;  // Stores weights for each principal component

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



Eigen::MatrixXd eigenVectorToMatrix(const Eigen::VectorXd& shapeVec) {
    Eigen::MatrixXd vertices(g_numVertices, 3);
    for (int v = 0; v < g_numVertices; ++v) {
        vertices.row(v) << shapeVec(3 * v), shapeVec(3 * v + 1), shapeVec(3 * v + 2);
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

    //// Apply deformation
    //Eigen::VectorXd deformation = g_eigenvectors.col(g_principalIndex) * g_weight;
    //Eigen::VectorXd deformedShape = baseShape + deformation;

    Eigen::VectorXd reconstructedShape = baseShape;

    // Apply all principal component weights
    for (int j = 0; j < g_eigenvectors.cols(); ++j) {
        reconstructedShape += g_weights[j] * g_eigenvectors.col(j);
    }

    // Update Polyscope mesh
    //std::vector<std::array<double, 3>> deformedPos = eigenVectorToVertices(deformedShape);
    std::vector<std::array<double, 3>> deformedPos = eigenVectorToVertices(reconstructedShape);
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




Eigen::MatrixXd generateAndVisualizePoints(const Eigen::MatrixXd& M, int numPoints=5000) {
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
    //int numPoints = 5000;

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
            //std::cout << "Generated Point: (" << point[0] << ", " << point[1] << ", " << point[2] << ")\n";

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

    //std::vector<int> selectedVerticesXXX = { 1 };
    //std::cout << "selectv:" << selectedVerticesXXX[0] << " " << selectedVerticesXXX.size() << std::endl;
    //std::vector<std::array<double, 3>> vertexColors1  = std::vector<std::array<double, 3>>(polyscopePoints.rows(), { {1.0, 1.0, 1.0} });
    //;
    //double radius = 0.005;
    //// Visualize the points in Polyscope
    ////auto obj1 = polyscope::registerPointCloud("Generated Points", polyscopePoints);
    //auto* obj1 = registerPointCloudWithColors("Point Cloud", polyscopePoints, radius, vertexColors1);
    //obj1->updatePointPositions(polyscopePoints); // Update Polyscope 
    //polyscope::state::userCallback = [&obj1, &polyscopePoints, &selectedVerticesXXX, &vertexColors1, &radius]() {

    //    std::cout << "selectv:" << selectedVerticesXXX[0] << " " << selectedVerticesXXX.size() << std::endl;

    //    HandlUserSelectionPCA(obj1, polyscopePoints, selectedVerticesXXX, vertexColors1, radius);
 
    //};
    //polyscope::show();


    return polyscopePoints;



}

Eigen::VectorXd extractInputPointsAsVector(const std::vector<int>& selectedVertices, const Eigen::MatrixXd& polyscopePoints) {
    Eigen::VectorXd shapeVec(3 * selectedVertices.size());

    for (size_t i = 0; i < selectedVertices.size(); ++i) {
        int v = selectedVertices[i];
        shapeVec.segment<3>(3 * i) = polyscopePoints.row(v).transpose();
    }

    return shapeVec;
}

Eigen::MatrixXd extractInputPointsAsMatrix(const std::vector<int>& selectedVertices, const Eigen::MatrixXd& polyscopePoints) {
    Eigen::MatrixXd pointsMatrix(selectedVertices.size(), 3);

    for (int i = 0; i < selectedVertices.size(); ++i) {
        int v = selectedVertices[i];

        // Access each component explicitly and assign to the matrix
        pointsMatrix(i, 0) = polyscopePoints(v, 0);  // x-coordinate
        pointsMatrix(i, 1) = polyscopePoints(v, 1);  // y-coordinate
        pointsMatrix(i, 2) = polyscopePoints(v, 2);  // z-coordinate


    }

    return pointsMatrix;
}


Eigen::VectorXd applyConstraints(
    const Eigen::VectorXd& meanShape,            // Mean shape (3N x 1)
    const Eigen::MatrixXd& principalComponents,  // Principal components (3N x k)
    const Eigen::VectorXi& constrainedIndices,   // Indices of constrained vertices (m x 1)
    const Eigen::MatrixXd& targetPositions)      // Target positions (m x 3)
{
    int numConstraints = constrainedIndices.size();
    int numVertices = meanShape.size() / 3;  // Number of vertices (N)

    // Ensure last column of PC is ignored, or leave it in depending on your needs
     //Eigen::MatrixXd PC = principalComponents.leftCols(principalComponents.cols() - 1); // Remove last PC
     Eigen::MatrixXd PC = principalComponents.leftCols(principalComponents.cols() - 3); // Remove last PC
    //Eigen::MatrixXd PC = principalComponents.col(0); // Retain all principal components

    // Selection matrix M (3m x 3N)
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(3 * numConstraints, 3 * numVertices);

    // Only modify the rows corresponding to constrained vertices in M
    for (int i = 0; i < numConstraints; ++i) {
        int idx = constrainedIndices(i);
        M.block(3 * i, 3 * idx, 3, 3) = Eigen::Matrix3d::Identity(); // Identity matrix for the constrained vertices
    }

    // Convert targetPositions to vector
    Eigen::VectorXd targetVec = Eigen::Map<const Eigen::VectorXd>(targetPositions.data(), 3 * numConstraints);

    // Compute the difference between target positions and the current positions of constrained points
    Eigen::VectorXd b = targetVec - M * meanShape;

    // Regularization: Add a small value to the diagonal to avoid instability
    double regularizationFactor = 1e-6;  // Regularization factor (tune this value)
    Eigen::MatrixXd regularizedMatrix = (M * PC).transpose() * (M * PC) + regularizationFactor * Eigen::MatrixXd::Identity((M * PC).cols(), (M * PC).cols());

    // Solve for the weights 'w' using a stable solver
    Eigen::VectorXd w = regularizedMatrix.ldlt().solve((M * PC).transpose() * b);

    // Compute the new shape by adding the weighted principal components
    Eigen::VectorXd newShape = meanShape + PC * w;

    // Returning the updated shape, keeping it flattened (3N x 1)
    return newShape;
}



Eigen::MatrixXd computeMeanShape(const std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>>& inputShapes) {
    if (inputShapes.empty()) {
        throw std::runtime_error("No shapes provided!");
    }

    int numShapes = inputShapes.size();
    int numPoints = inputShapes[0].first.rows();  // Number of points in each shape
    int dim = inputShapes[0].first.cols();       // Should be 3 (x, y, z)

    // Initialize the matrix to accumulate the sum of all shapes
    Eigen::MatrixXd sumShape = Eigen::MatrixXd::Zero(numPoints, dim);

    // Sum all shape matrices
    for (const auto& shape : inputShapes) {
        sumShape += shape.first;  // Add the vertex positions (x, y, z)
    }

    // Divide by the number of shapes to compute the mean
    Eigen::MatrixXd meanShape = sumShape / numShapes;

    return meanShape;
}

Eigen::VectorXd convert2DOffsetsTo3D(const Eigen::VectorXd& cFFD_2D)
{
    // The input cFFD_2D has size 2N
    int N = static_cast<int>(cFFD_2D.size() / 2);
    // We'll create a 3N vector for 3D offsets
    Eigen::VectorXd cFFD_3D(3 * N);

    for (int i = 0; i < N; ++i)
    {
        double offsetX = cFFD_2D(2 * i + 0);
        double offsetY = cFFD_2D(2 * i + 1);

        cFFD_3D(3 * i + 0) = offsetX;
        cFFD_3D(3 * i + 1) = offsetY;
        cFFD_3D(3 * i + 2) = 0.0;    // no offset in z
    }

    return cFFD_3D;
}


#include <Eigen/Core>
#include <Eigen/Dense>
#include <tuple>
#include <vector>
#include <iostream>
#include <limits>
#include <cmath>
#include <algorithm>

// A small struct storing how a single vertex is influenced by 4 control points:
struct CPInfluence
{
    // Indices of the 4 corner control points (flattened)
    int cp00, cp10, cp01, cp11;
    // Bilinear weights for each corner
    double w00, w10, w01, w11;
};

// Compute a bounding box that includes both the shape (V) and the selected 2D points.
static void unifiedBoundingBox(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXd& selected_points,
    double& min_x, double& max_x,
    double& min_y, double& max_y)
{
    min_x = std::min(V.col(0).minCoeff(), selected_points.col(0).minCoeff());
    max_x = std::max(V.col(0).maxCoeff(), selected_points.col(0).maxCoeff());
    min_y = std::min(V.col(1).minCoeff(), selected_points.col(1).minCoeff());
    max_y = std::max(V.col(1).maxCoeff(), selected_points.col(1).maxCoeff());
    // If the shape and points are identical or extremely close, ensure non-zero range:
    if (max_x - min_x < 1e-12) { max_x = min_x + 1.0; }
    if (max_y - min_y < 1e-12) { max_y = min_y + 1.0; }
}

// Precompute how each vertex is influenced by the grid control points.
// For each vertex i, we find which cell of the grid it falls into, compute
// bilinear weights (w00, w10, w01, w11), and store the 4 corner CP indices.
static std::vector<CPInfluence> precomputeFFDWeights(
    const Eigen::MatrixXd& V,
    double min_x, double max_x,
    double min_y, double max_y,
    int gridX, int gridY)
{
    const int N = (int)V.rows();
    std::vector<CPInfluence> influences(N);

    for (int i = 0; i < N; i++)
    {
        double x = V(i, 0);
        double y = V(i, 1);

        // Map (x,y) into grid space in [0, gridX-1], [0, gridY-1]
        double sx = (x - min_x) / (max_x - min_x) * (gridX - 1);
        double sy = (y - min_y) / (max_y - min_y) * (gridY - 1);

        int ix = (int)std::floor(sx);
        int iy = (int)std::floor(sy);

        // Clamp to avoid out-of-bounds
        if (ix < 0) ix = 0;
        if (iy < 0) iy = 0;
        if (ix >= gridX - 1) ix = gridX - 2;
        if (iy >= gridY - 1) iy = gridY - 2;

        double u = sx - ix;
        double v = sy - iy;

        // Flattened CP indices
        int cp00 = ix * gridY + iy;
        int cp10 = (ix + 1) * gridY + iy;
        int cp01 = ix * gridY + (iy + 1);
        int cp11 = (ix + 1) * gridY + (iy + 1);

        CPInfluence inf;
        inf.cp00 = cp00; inf.cp10 = cp10; inf.cp01 = cp01; inf.cp11 = cp11;
        inf.w00 = (1 - u) * (1 - v);
        inf.w10 = u * (1 - v);
        inf.w01 = (1 - u) * v;
        inf.w11 = u * v;
        influences[i] = inf;
    }
    return influences;
}

// Deform the entire shape (in x,y) using the precomputed FFD weights and the param vector.
// param.size() = 2*(gridX*gridY), storing (dx, dy) for each CP in row-major flattening.
static Eigen::MatrixXd deformShape(
    const Eigen::MatrixXd& V,
    const std::vector<CPInfluence>& influences,
    const Eigen::VectorXd& param)
{
    Eigen::MatrixXd newV = V;
    for (int i = 0; i < (int)V.rows(); i++)
    {
        const auto& inf = influences[i];

        // Retrieve the displacements from param
        auto getDisp = [&](int cpIndex) {
            double dx = param(2 * cpIndex + 0);
            double dy = param(2 * cpIndex + 1);
            return std::make_pair(dx, dy);
        };

        auto [dx00, dy00] = getDisp(inf.cp00);
        auto [dx10, dy10] = getDisp(inf.cp10);
        auto [dx01, dy01] = getDisp(inf.cp01);
        auto [dx11, dy11] = getDisp(inf.cp11);

        // Weighted sum
        double dx = inf.w00 * dx00 + inf.w10 * dx10 + inf.w01 * dx01 + inf.w11 * dx11;
        double dy = inf.w00 * dy00 + inf.w10 * dy10 + inf.w01 * dy01 + inf.w11 * dy11;

        newV(i, 0) = V(i, 0) + dx;
        newV(i, 1) = V(i, 1) + dy;
        // z stays the same
    }
    return newV;
}

// The main function performing "non-rigid ICP" with a linear solve each iteration.
std::tuple<Eigen::MatrixXd, Eigen::VectorXd, std::vector<int>>
nonRigidICPFFD_Linear(
    const Eigen::MatrixXd& selected_points,  // Mx2
    const Eigen::MatrixXd& V,               // Nx3
    int gridX = 10, int gridY = 10,
    int outerIterations = 100)
{
    // 1) Compute a bounding box that covers both shape + selected points
    double min_x, max_x, min_y, max_y;
    unifiedBoundingBox(V, selected_points, min_x, max_x, min_y, max_y);

    // 2) Precompute how each vertex is influenced by the grid CPs
    std::vector<CPInfluence> influences =
        precomputeFFDWeights(V, min_x, max_x, min_y, max_y, gridX, gridY);

    // 3) Initialize param = zero displacement
    int numCP = gridX * gridY;
    Eigen::VectorXd param = Eigen::VectorXd::Zero(2 * numCP);

    // We'll keep track of correspondences: for each 2D point, which vertex?
    int M = (int)selected_points.rows();
    std::vector<int> correspondences(M, 0);

    // We also maintain a "deformed shape" for each iteration
    Eigen::MatrixXd deformedV = V;

    // 4) Outer iteration loop
    for (int iter = 0; iter < outerIterations; iter++)
    {
        // (A) Deform the shape with the current param
        deformedV = deformShape(V, influences, param);

        // (B) Reassign correspondences: each 2D point -> nearest DEFORMED vertex in x,y
        for (int m = 0; m < M; m++)
        {
            double px = selected_points(m, 0);
            double py = selected_points(m, 1);

            double bestDist = std::numeric_limits<double>::max();
            int bestIdx = -1;
            for (int i = 0; i < (int)deformedV.rows(); i++)
            {
                double dx = deformedV(i, 0) - px;
                double dy = deformedV(i, 1) - py;
                double dist2 = dx * dx + dy * dy;
                if (dist2 < bestDist)
                {
                    bestDist = dist2;
                    bestIdx = i;
                }
            }
            correspondences[m] = bestIdx;
        }

        // (C) Build the linear system A * param = b, with 2*M equations
        //     Because for each match:
        //       x_i + sum_c (w_i,c * dx_c) = p_x
        //       y_i + sum_c (w_i,c * dy_c) = p_y
        //     We'll do row 2*m for the x eq, row 2*m+1 for y eq.
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2 * M, 2 * numCP);
        Eigen::VectorXd b = Eigen::VectorXd::Zero(2 * M);

        for (int m = 0; m < M; m++)
        {
            int vidx = correspondences[m];
            const auto& inf = influences[vidx];

            double vx = V(vidx, 0);
            double vy = V(vidx, 1);
            double px = selected_points(m, 0);
            double py = selected_points(m, 1);

            // For x eq: row = 2*m
            double rx = px - vx; // we want sum_c w_i,c * dx_c = rx
            // Fill columns for the 4 corners
            auto setX = [&](int cpIndex, double w) {
                A(2 * m + 0, 2 * cpIndex + 0) += w; // dx component
            };

            setX(inf.cp00, inf.w00);
            setX(inf.cp10, inf.w10);
            setX(inf.cp01, inf.w01);
            setX(inf.cp11, inf.w11);

            b(2 * m + 0) = rx;

            // For y eq: row = 2*m+1
            double ry = py - vy;
            auto setY = [&](int cpIndex, double w) {
                A(2 * m + 1, 2 * cpIndex + 1) += w; // dy component
            };

            setY(inf.cp00, inf.w00);
            setY(inf.cp10, inf.w10);
            setY(inf.cp01, inf.w01);
            setY(inf.cp11, inf.w11);

            b(2 * m + 1) = ry;
        }

        // Solve the linear system for param (in least-squares sense)
        Eigen::VectorXd newParam = A.colPivHouseholderQr().solve(b);

        // (Optional) Check how much param changes, or we can just assign
        param = newParam;
    }

    // One last deformation with the final param
    deformedV = deformShape(V, influences, param);

    return std::make_tuple(deformedV, param, correspondences);
}







// Function to compute optimal weights w using gradient descent with Jacobian
std::pair<Eigen::VectorXd, std::vector<int>> computeOptimalWeightsWithGD(
    const Eigen::VectorXd& g_meanShape,      // (3N, 1) mean shape stored in a flat format
    const Eigen::MatrixXd& g_eigenvectors,  // (3N, k) eigenvectors stored in a flat format
    const Eigen::MatrixXd& selected_points, // (M, 2) 2D outline points stored as [(x0,y0), (x1,y1), ...]
    double learning_rate = 0.01,          // Learning rate for gradient descent
    double tolerance = 1e-6               // Convergence threshold

) {
    int N = g_meanShape.size() / 3;  // Number of 3D points
    int k = g_eigenvectors.cols();   // Number of principal components
    int M = selected_points.rows();  // Number of outline points

    // Step 1: Extract x, y coordinates from g_meanShape (size 2N × 1)
    Eigen::VectorXd proj_meanShape(2 * N);
    for (int i = 0; i < N; ++i) {
        proj_meanShape(2 * i) = g_meanShape(3 * i);         // x_i
        proj_meanShape(2 * i + 1) = g_meanShape(3 * i + 1); // y_i
    }

    // Step 2: Extract x, y components from eigenvectors (size 2N × k)
    Eigen::MatrixXd proj_eigenVectors(2 * N, k);
    for (int i = 0; i < N; ++i) {
        proj_eigenVectors.row(2 * i) = g_eigenvectors.row(3 * i);         // x components
        proj_eigenVectors.row(2 * i + 1) = g_eigenvectors.row(3 * i + 1); // y components
    }

    // Step 3: Flatten selected_points (M, 2) ? (2M, 1)
    Eigen::VectorXd selected_points_flat(2 * M);
    for (int j = 0; j < M; ++j) {
        selected_points_flat(2 * j) = selected_points(j, 0); // x_j
        selected_points_flat(2 * j + 1) = selected_points(j, 1); // y_j
    }

    // Step 4: Find closest mean shape point for each outline point
    std::vector<int> closestIndices(M);
    for (int j = 0; j < M; ++j) {
        Eigen::Vector2d outline_pt = selected_points.row(j);
        double minDist = std::numeric_limits<double>::max();
        int best_idx = -1;

        for (int i = 0; i < N; ++i) {
            Eigen::Vector2d shape_pt = proj_meanShape.segment<2>(2 * i);
            double dist = (shape_pt - outline_pt).squaredNorm();
            if (dist < minDist) {
                minDist = dist;
                best_idx = i;
            }
        }
        closestIndices[j] = best_idx;
    }

    // Step 5: Gradient Descent for optimal weights w using Jacobian
    Eigen::VectorXd w = Eigen::VectorXd::Zero(k);  // Initialize weights (k × 1)
    double prev_energy = std::numeric_limits<double>::max(); // Previous iteration's energy
    int iteration = 0;
    Eigen::VectorXd reconstructedShape;
    // Step 6: Initialize reconstructed shape with mean shape
    Eigen::VectorXd proj_reconstructedShape = proj_meanShape;

    for (int iteration = 0; iteration < 100000; ++iteration) {
        // Step 5.1: Compute residuals using the latest projected reconstructed shape
        Eigen::VectorXd residuals(2 * M);
        for (int j = 0; j < M; ++j) {
            int idx = closestIndices[j];
            residuals.segment<2>(2 * j) = selected_points_flat.segment<2>(2 * j) - proj_reconstructedShape.segment<2>(2 * idx);
        }

        // Step 5.2: Compute Jacobian
        Eigen::MatrixXd jacobian(2 * M, k);
        for (int j = 0; j < M; ++j) {
            int idx = closestIndices[j];
            jacobian.block(2 * j, 0, 2, k) = proj_eigenVectors.block(2 * idx, 0, 2, k);
        }

        // Step 5.3: Compute gradient and Hessian
        Eigen::VectorXd gradient = jacobian.transpose() * residuals;
        Eigen::MatrixXd hessian = jacobian.transpose() * jacobian;

        // Step 5.4: Compute update step
        Eigen::VectorXd update = hessian.ldlt().solve(gradient);
        w += learning_rate * update;

        // Step 6: Update the reconstructed shape with the new weights
        proj_reconstructedShape = proj_meanShape;  // Reset to mean shape
        for (int j = 0; j < k; ++j) {
            proj_reconstructedShape += w[j] * proj_eigenVectors.col(j);
        }

        // Step 7: Compute energy
        double energy = residuals.squaredNorm();
        std::cout << "Iteration " << iteration << ", Energy: " << energy << std::endl;

        // Stopping condition: Check if the update is small
        if (update.norm() < 1e-6) {
            std::cout << "Converged at iteration " << iteration << std::endl;
            break;
        }
    }

    return std::make_pair(w, closestIndices); // Return the optimized weights
}


// Function for Laplacian smoothing
void laplacianSmoothing(Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int iterations = 10) {
    int numVertices = V.rows();
    int numFaces = F.rows();

    // Create a matrix to store vertex neighbors
    std::vector<std::vector<int>> vertexNeighbors(numVertices);

    // Step 1: Build the adjacency list (neighbors of each vertex)
    for (int f = 0; f < numFaces; ++f) {
        for (int i = 0; i < 3; ++i) {
            int v1 = F(f, i);
            int v2 = F(f, (i + 1) % 3);
            vertexNeighbors[v1].push_back(v2);
            vertexNeighbors[v2].push_back(v1);
        }
    }

    // Step 2: Perform Laplacian smoothing
    for (int it = 0; it < iterations; ++it) {
        Eigen::MatrixXd V_new = V;  // Copy of the current vertices

        // Iterate through each vertex and update its position
        for (int v = 0; v < numVertices; ++v) {
            Eigen::RowVector3d newPos(0, 0, 0);
            int count = 0;

            // Average the position of all neighbors
            for (int neighbor : vertexNeighbors[v]) {
                newPos += V.row(neighbor);
                count++;
            }

            // Set the new vertex position as the average
            if (count > 0) {
                newPos /= count;
                V_new.row(v) = newPos;
            }
        }

        // Update the vertex positions with the smoothed version
        V = V_new;
    }
}


#include <igl/arap.h>
// Function to perform ARAP on the reconstructed shape
Eigen::MatrixXd performARAP(
    const Eigen::MatrixXd& V_init,   // (N, 3) Initial vertices (from PCA)
    const Eigen::MatrixXi& F,        // (M, 3) Faces of the mesh
    const Eigen::MatrixXd& V_fixed,  // (M, 2) 2D outline points (used for constraints)
    const std::vector<int>& fixedIndices // Indices of points to keep fixed
) {

    std::cout << "ARAP" << std::endl;
    // Print input sizes for debugging
    std::cout << "V_init size: " << V_init.rows() << " x " << V_init.cols() << std::endl;
    std::cout << "F size: " << F.rows() << " x " << F.cols() << std::endl;
    std::cout << "V_fixed size: " << V_fixed.rows() << " x " << V_fixed.cols() << std::endl;
    std::cout << "fixedIndices size: " << fixedIndices.size() << std::endl;
    int dim = 3; // ARAP in 3D

    // Step 1: Convert std::vector<int> fixedIndices to Eigen::MatrixXi
    Eigen::MatrixXi fixedIndicesMat(fixedIndices.size(), 1);
    for (int i = 0; i < fixedIndices.size(); ++i) {
        fixedIndicesMat(i, 0) = fixedIndices[i];
    }

    // Step 2: Setup ARAP energy
    igl::ARAPData arap_data;
    arap_data.with_dynamics = false; // Standard ARAP
    igl::arap_precomputation(V_init, F, dim, fixedIndicesMat, arap_data);

    // Step 3: Initialize ARAP solver with initial shape
    Eigen::MatrixXd V = V_init;

    // Step 4: Set constraint points
    Eigen::MatrixXd bc(fixedIndices.size(), 3);
    for (int i = 0; i < fixedIndices.size(); ++i) {
        bc.row(i).head(2) = V_fixed.row(i); // Keep x, y from selected outline
        bc(i, 2) = V(fixedIndices[i], 2);   // Let ARAP solve for Z
    }

    // Step 5: Solve iterative ARAP
    igl::arap_solve(bc, arap_data, V);

    return V; // Final deformed shape
}

void printThicknesses(const std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>>& inputShapes) {
    for (size_t i = 0; i < inputShapes.size(); ++i) {
        const Eigen::MatrixXd& V = inputShapes[i].first;

        if (V.rows() == 0 || V.cols() < 3) {
            std::cout << "Shape " << i << " has invalid vertex data.\n";
            continue;
        }

        double minZ = V.col(2).minCoeff();
        double maxZ = V.col(2).maxCoeff();
        double thickness = std::abs(maxZ - minZ);

        std::cout << "Shape " << i << " thickness (Z): " << thickness << std::endl;
    }
}

double computeThickness(const Eigen::MatrixXd& V) {
    if (V.rows() == 0 || V.cols() < 3) {
        std::cerr << "Invalid vertex data.\n";
        return 0.0;
    }

    double minZ = V.col(2).minCoeff();
    double maxZ = V.col(2).maxCoeff();
    auto val = std::abs(maxZ - minZ);
    std::cout << "mean shape thickness: " << val << std::endl;

    return val;
}



// Main PCA computation and visualization setup
void performPCAAndEditWithVisualization(const std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>>& inputShapes) {
    if (inputShapes.empty()) {
        std::cerr << "No shapes provided." << std::endl;
        return;
    }

    printThicknesses(inputShapes);

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
    g_weights.resize(g_eigenvectors.cols(), 0.0f); // Initialize weights to zero


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


    Eigen::MatrixXd meanShape3dMat = computeMeanShape(inputShapes);
    computeThickness(meanShape3dMat);
    //Eigen::MatrixXd polyscopePoints = meanShape3dMat
    int defaultGridTotal = 10000;
    Eigen::MatrixXd polyscopePoints = generateAndVisualizePoints(inputShapes[0].first, defaultGridTotal);

    //std::vector<int> selectedVerticesXXX = { 1 };
    //std::vector<int> selectedVerticesXXX = {1634, 2430, 3315};
    //std::vector<int> selectedVerticesXXX = { 2768, 2769, 2767, 2972, 3041, 3110, 2758, 2689, 2690, 2624, 2625, 2558, 2559, 2491, 2420, 2419, 2418, 2347, 2346, 2635, 2636, 2637, 2639, 2570, 2501, 2431, 2362, 2290, 2219, 2218, 2216, 2210, 2213, 2214, 2215, 1921, 1850, 2344, 2209, 2208, 2207, 2136, 2064, 2626, 2622, 3245, 3043, 2904, 2976, 2978, 2837, 2908, 3246, 3247, 3248, 3179, 3036, 3105, 3174, 2966, 2897, 2828, 2273, 2202, 2059, 1987, 1916, 1845, 1779, 1708, 1774, 1703, 1632, 1637, 1566, 1495, 1353, 1424, 1351, 1282, 1420, 1561, 1490, 1992, 2060, 2131 };
    //std::vector<int> selectedVerticesXXX = { 3835, 3936, 4037, 4138, 4239, 4340, 4341, 4442, 4443, 4444, 4546, 4445, 4648, 4750, 4751, 4649, 4852, 4853, 4953, 4954, 4955, 4856, 4959, 5060, 4958, 5159, 5158, 5157, 5255, 5254, 5252, 5349, 5350, 5446, 5347, 5544, 5643, 5642, 5641, 5640, 5639, 5738, 5838, 5938, 6038, 6138, 6237, 6337, 6437, 6537, 6638, 6639, 6640, 6642, 6543, 6444, 6345, 6244, 6145, 6045, 5946, 5947, 5848, 5750, 5849, 5652, 5653, 5654, 5655, 5656, 5657, 5658, 5659, 5660, 5561, 5462, 5363, 5364, 5365, 5266, 5267, 5268, 5270, 5269, 5171, 5071, 4971, 4970, 4868, 4867, 4766, 4765, 4663, 4662, 4661, 4660, 4558, 4557, 4556, 4454, 4453, 4352, 4351, 4249, 4148, 4147, 4046, 3844, 3945, 3743, 3641, 3642, 3539, 3438, 3437, 3536, 3635, 3735, 4555 }; // 2 blader
    std::vector<int> selectedVerticesXXX = {3937, 4038, 4139, 4240, 4341, 4442, 4443, 4544, 4545, 4646, 4647, 4748, 4749, 4750, 4751, 4752, 4753, 4754, 4755, 4756, 4757, 4858, 4958, 5058, 5157, 5156, 5155, 5254, 5253, 5252, 5251, 5350, 5349, 5448, 5447, 5546, 5545, 5644, 5743, 5842, 5941, 6040, 6139, 6238, 6337, 6436, 6536, 6636, 6736, 6835, 6934, 7033, 7132, 7232, 7333, 7434, 7435, 7436, 7337, 7338, 7239, 7140, 7041, 6942, 6843, 6743, 6644, 6544, 6445, 6345, 6245, 6146, 6046, 5947, 5848, 5749, 5651, 5652, 5653, 5750, 5654, 5655, 5656, 5657, 5558, 5559, 5561, 5560, 5462, 5463, 5364, 5265, 5166, 5066, 4966, 4866, 4766, 4665, 4564, 4463, 4362, 4361, 4360, 4359, 4358, 4357, 4256, 4255, 4254, 4253, 4252, 4251, 4250, 4249, 4248, 4147, 4046, 3945, 3844, 3744, 3643, 3543, 3442, 3342, 3242, 3142, 3042, 2942, 2841, 2840, 2839, 2738, 2737, 2736, 2735, 2733, 2832, 2734, 2931, 3031, 3131, 3232, 3333, 3434, 3535, 3635, 3736, 3836, 3837, 4957, 5057, 4857 }; // 2 blader
    //std::vector<int> selectedVerticesXXX = { 4146, 4347, 4447, 4548, 4649, 4750, 4851, 4852, 4853, 4854, 4855, 4856, 4757, 4758, 4660, 4562, 4563, 4465, 4464, 4466, 4567, 4667, 4768, 4868, 5068, 5268, 5267, 5266, 5165, 5164, 5163, 5062, 5061, 5160, 5159, 5257, 5158, 5455, 5356, 5454, 5352, 5351, 5249, 5248, 5247, 5246, 5344, 5443, 5542, 5741, 5942, 5841, 6142, 6042, 6241, 6441, 6341, 6541, 6739, 6640, 6838, 6837, 6836, 6835, 6734, 6534, 6434, 6335, 6236, 6136, 5937, 5837, 6036, 5736, 5636, 5437, 5536, 5338, 5240, 5141, 5042, 4943, 4843, 4743, 4643, 4542, 4442, 4341, 4240, 4239, 4138, 4037, 4036, 3935, 3934, 3933, 3932, 3931, 3930, 3730, 3630, 3531, 3432, 3333, 3234, 3235, 3236, 3337, 3438, 3538, 3739, 3638, 3840, 3941, 3942, 4044, 4045, 6334, 6634 }; // 3blader

    std::cout << "selectv:" << selectedVerticesXXX[0] << " " << selectedVerticesXXX.size() << std::endl;
    std::vector<std::array<double, 3>> vertexColors1  = std::vector<std::array<double, 3>>(polyscopePoints.rows(), { {1.0, 1.0, 1.0} });

    double radius = 0.005;
    auto* obj1 = registerPointCloudWithColors("Point Cloud", polyscopePoints, radius, vertexColors1);
    obj1->updatePointPositions(polyscopePoints); // Update Polyscope 
    //

    // Register the mean shape
    polyscope::registerSurfaceMesh("Mean Shape", eigenVectorToVertices(g_meanShape), g_faceList);
    polyscope::registerSurfaceMesh("Deformed Shape", eigenVectorToVertices(g_meanShape), g_faceList);
    g_deformedMesh = polyscope::getSurfaceMesh("Deformed Shape");

    
    // Global or static storage for optimal weights
    static Eigen::VectorXd optimalWeights;
    static Eigen::VectorXd def;
    std::vector<int> closestIndices;
    static float lambda = 0.1f; // Default regularization strength

    // User controls via ImGui
    polyscope::state::userCallback = [&obj1, &polyscopePoints, &selectedVerticesXXX, &vertexColors1, &radius, &inputShapes, &closestIndices]() {
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

        // Sliders for each principal component weight
        for (int i = 0; i < totalModes; ++i) {
            ImGui::SliderFloat(("Weight " + std::to_string(i)).c_str(), &g_weights[i], -10000.0f, 10000.0f);
        }

        // Reset all weights button
        if (ImGui::Button("Reset All Weights")) {
            std::fill(g_weights.begin(), g_weights.end(), 0.0f);
        }

        // Generate random sample button
        if (ImGui::Button("Generate Random Sample")) {
            generateRandomShape();
        }

        // Update deformed shape
        updateDeformedMesh();

        // Print selected points
        if (ImGui::Button("Print Selected Points")) {
            std::cout << "Selected Points Coordinates:\n";
            for (int idx : selectedVerticesXXX) {
                if (idx >= 0 && idx < polyscopePoints.rows()) {
                    std::cout << "Index " << idx << ": (" << polyscopePoints(idx, 0) << ", " << polyscopePoints(idx, 1) << ", " << polyscopePoints(idx, 2) << ")\n";
                }
                else {
                    std::cout << "Invalid index: " << idx << "\n";
                }
            }
        }
        if (ImGui::Button("Clear Selected Vertices")) {
            selectedVerticesXXX.clear();
            std::cout << "Selected vertices list cleared." << std::endl;
        }

        if (ImGui::Button("Show Selected Points (Green)")) {
            Eigen::MatrixXd selectedPointsGreen(selectedVerticesXXX.size(), 3);
            int validCount = 0;

            for (int idx : selectedVerticesXXX) {
                if (idx >= 0 && idx < polyscopePoints.rows()) {
                    selectedPointsGreen.row(validCount) = polyscopePoints.row(idx);
                    ++validCount;
                }
                else {
                    std::cout << "Invalid index: " << idx << "\n";
                }
            }

            selectedPointsGreen.conservativeResize(validCount, Eigen::NoChange);

            if (validCount > 0) {
                std::vector<std::array<double, 3>> greenColors(validCount, { {1.0, 0.0, 0.0} });
                registerPointCloudWithColors("Selected Points (Red)", selectedPointsGreen, 0.01, greenColors);
            }
        }


        HandlUserSelectionPCA(obj1, polyscopePoints, selectedVerticesXXX, vertexColors1, radius);
        
        // Button to compute optimal weights
        if (ImGui::Button("Compute Optimal Weights")) {
            std::cout << "Computing Optimal Weights..." << std::endl;
            Eigen::MatrixXd selected_points = extractInputPointsAsMatrix(selectedVerticesXXX, polyscopePoints);
            auto result = computeOptimalWeightsWithGD(g_meanShape, g_eigenvectors, selected_points);
            optimalWeights = result.first;
            closestIndices = result.second;

            std::cout << "Optimal Weights Computed:\n" << optimalWeights << std::endl;
        }




        // Slider for custom lambda
        ImGui::SliderFloat("Lambda", &lambda, 0.0f, 1.0f, "%.01f");

        // Button to apply stored weights with lambda regularization
        if (ImGui::Button("Apply Weights with Custom Lambda")) {
            std::cout << "Applying stored optimal weights with lambda = " << lambda << std::endl;

            Eigen::VectorXd reconstructedShape = g_meanShape;

            for (int j = 0; j < g_eigenvectors.cols() - 1; ++j) {
                reconstructedShape += optimalWeights[j]  * lambda * g_eigenvectors.col(j);
            }



            polyscope::registerSurfaceMesh("Fitted shape", eigenVectorToVertices(reconstructedShape), g_faceList);


            //Eigen::VectorXd reconstructedShape2 = reconstructedShape;
            //Eigen::VectorXd cFFD_3D = convert2DOffsetsTo3D(def); // size 3N
            //reconstructedShape2 += cFFD_3D;




            //Eigen::MatrixXd V_arap = performARAP(eigenVectorToMatrix(reconstructedShape), inputShapes[0].second, extractInputPointsAsMatrix(selectedVerticesXXX, polyscopePoints), closestIndices);

            //polyscope::registerSurfaceMesh("Fitted shape ARAP", V_arap, inputShapes[0].second);


            //Eigen::MatrixXd V_lap = V_arap;
            //Eigen::MatrixXi F_lap = inputShapes[0].second;
            //laplacianSmoothing(V_lap, F_lap);

            //polyscope::registerSurfaceMesh("Fitted shape ARAP + Laplacian", V_lap, F_lap);

        }

        // Button to overwrite polyscopePoints

        static int InputInt = 5000; // Default value for input
        ImGui::SliderInt("Input Points", &InputInt, 1000, 50000); // Adjust range as needed

        // Snap InputInt to the nearest 100
        InputInt = (InputInt / 100) * 100;

        if (ImGui::Button("Generate New Points Grid")) {
            std::cout << "Generating new polyscope points with InputInt = " << InputInt << std::endl;

            polyscopePoints = generateAndVisualizePoints(inputShapes[0].first, InputInt);

            // Clear selection to avoid invalid indices
            selectedVerticesXXX.clear();
            std::cout << "Cleared selected vertices due to grid change." << std::endl;

            // Re-register point cloud to avoid inconsistencies
            polyscope::removePointCloud("Point Cloud");
            registerPointCloudWithColors("Point Cloud", polyscopePoints, radius, vertexColors1);
        }


        // Grid size sliders
        static int GridX = 10; // Default grid size X
        static int GridY = 10; // Default grid size Y

        ImGui::SliderInt("GridX", &GridX, 2, 50); // Adjust range as needed
        ImGui::SliderInt("GridY", &GridY, 2, 50); // Adjust range as needed




        if (ImGui::Button("Compute FFD")) {
            std::cout << "Computing FFD" << std::endl;
            Eigen::MatrixXd selected_points = extractInputPointsAsMatrix(selectedVerticesXXX, polyscopePoints);

            Eigen::VectorXd reconstructedShape = g_meanShape;

            for (int j = 0; j < g_eigenvectors.cols() - 1; ++j) {
                reconstructedShape += optimalWeights[j] * lambda * g_eigenvectors.col(j);
            }

            auto res2 = nonRigidICPFFD_Linear(selected_points, eigenVectorToMatrix(reconstructedShape), GridX, GridY);

            Eigen::MatrixXd V_new = std::get<0>(res2);

            polyscope::registerSurfaceMesh("A", eigenVectorToMatrix(reconstructedShape), inputShapes[0].second);
            polyscope::registerSurfaceMesh("FFD result", V_new, inputShapes[0].second);

        }
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





