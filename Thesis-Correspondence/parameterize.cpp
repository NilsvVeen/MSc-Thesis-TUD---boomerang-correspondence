
#include "parameterize.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include "stl_utils.h"



// Function to reverse the order of rows in an Eigen::MatrixXd
Eigen::MatrixXd reverseOrder(const Eigen::MatrixXd& V) {
    // Create a new matrix to hold the reversed rows
    Eigen::MatrixXd reversed(V.rows(), V.cols());

    for (int i = 0; i < V.rows(); ++i) {
        reversed.row(i) = V.row(V.rows() - 1 - i); // Assign the row from the end to the new matrix
    }

    return reversed;
}

// Function to compute Euclidean distance between two points
double computeDistance(const Eigen::VectorXd& p1, const Eigen::VectorXd& p2) {
    return (p1 - p2).norm(); // Calculate the Euclidean distance
}

// Function to compute the centroid of the vertices
Eigen::Vector2d computeCentroid(const Eigen::MatrixXd& V) {
    return V.leftCols<2>().colwise().mean(); // Mean of x and y coordinates
}

// Function to sort vertices based on proximity and ensure clockwise ordering in 2D (x, y)
Eigen::MatrixXd sortVerticesByProximity(const Eigen::MatrixXd& V) {
    int numVertices = V.rows();
    Eigen::MatrixXd sortedVertices(numVertices, V.cols()); // To store the sorted vertices

    // Compute the centroid of the vertices (average position)
    Eigen::Vector2d centroid = computeCentroid(V);

    // Lambda function to compute the angle of a vertex relative to the centroid
    auto computeAngle = [&](const Eigen::Vector2d& vertex) -> double {
        Eigen::Vector2d relativePos = vertex - centroid;
        return std::atan2(relativePos.y(), relativePos.x()); // Angle in radians
    };

    // Create a vector of pairs (angle, vertex index) for sorting
    std::vector<std::pair<double, int>> angleIndexPairs;
    for (int i = 0; i < numVertices; ++i) {
        Eigen::Vector2d vertex2D = V.row(i).head<2>(); // Get x, y coordinates
        double angle = computeAngle(vertex2D);
        angleIndexPairs.emplace_back(angle, i); // Store the angle and corresponding index
    }

    // Sort the vertices based on the angle (in ascending order for clockwise)
    std::sort(angleIndexPairs.begin(), angleIndexPairs.end(), [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
        return a.first < b.first; // Sort by angle (clockwise)
        });

    // Use the sorted indices to reorder the vertices
    for (int i = 0; i < numVertices; ++i) {
        sortedVertices.row(i) = V.row(angleIndexPairs[i].second); // Add the sorted vertex
    }

    return sortedVertices;
}

// Function to compute the Euclidean distance between two points (x,y)
//double distance2d(const Point3& a, const Point3& b) {
//    double dx = b.x - a.x;
//    double dy = b.y - a.y;
//    return std::sqrt(dx * dx + dy * dy);
//}

double distance2d(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2) {
    return (v1.head<2>() - v2.head<2>()).norm(); // Using the first two components (x, y)
}


std::vector<double> computeUnitParametrization(const Eigen::MatrixXd& vertices) {
    int n = vertices.rows();
    if (n < 2) {
        throw std::invalid_argument("At least two vertices are required to form a curve.");
    }

    std::vector<double> params(n, 0.0);
    double lengthCurve = 0.0;
    std::vector<double> segments(n - 1, 0.0);

    for (int i = 0; i < n - 1; ++i) {
        segments[i] = distance2d(vertices.row(i), vertices.row(i + 1));
        lengthCurve += segments[i];
    }

    // Compute parameter values for each vertex
    for (int i = 1; i < n; ++i) {
        params[i] = params[i - 1] + (segments[i - 1] / lengthCurve);
    }

    // Ensure the last vertex is mapped to 1
    params[n - 1] = 1.0;

    return params;
}





// Downsample points 

// Function to downsample a point cloud to match the number of points in the smaller point cloud 
// note could be more fancy
Eigen::MatrixXd downsamplePointCloud(const Eigen::MatrixXd& V, int targetNumPoints) {
    int originalNumPoints = V.rows();

    // If the number of points is already smaller or equal to the target, return the original
    if (originalNumPoints <= targetNumPoints) {
        return V;
    }

    // Create a new matrix to store the downsampled points
    Eigen::MatrixXd downsampledV(targetNumPoints, V.cols());

    // Calculate step size for sampling
    double step = static_cast<double>(originalNumPoints - 1) / (targetNumPoints - 1);

    // Sample points at regular intervals
    for (int i = 0; i < targetNumPoints; ++i) {
        int index = static_cast<int>(i * step);
        downsampledV.row(i) = V.row(index);
    }

    return downsampledV;
}

// Function to equalize the number of points in both point clouds
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> equalizePointClouds(const Eigen::MatrixXd& V1, const Eigen::MatrixXd& V2) {
    int numPoints1 = V1.rows();
    int numPoints2 = V2.rows();

    // Determine the minimum number of points
    int targetNumPoints = std::min(numPoints1, numPoints2);

    // Downsample both point clouds to have the same number of points
    Eigen::MatrixXd V1_downsampled = downsamplePointCloud(V1, targetNumPoints);
    Eigen::MatrixXd V2_downsampled = downsamplePointCloud(V2, targetNumPoints);

    return std::make_pair(V1_downsampled, V2_downsampled);
}

// Main function to downsample and return equalized sorted vertices
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> getEqualizedPointClouds(const Eigen::MatrixXd& sortedVertices1, const Eigen::MatrixXd& sortedVertices2) {

    // Equalize the number of points in both point clouds
    auto [equalizedV1, equalizedV2] = equalizePointClouds(sortedVertices1, sortedVertices2);

    // Return the downsampled point clouds
    return { equalizedV1, equalizedV2 };
}
