
#include "parameterize.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include "stl_utils.h"



// Function to compute Euclidean distance between two points
double computeDistance(const Eigen::VectorXd& p1, const Eigen::VectorXd& p2) {
    return (p1 - p2).norm(); // Calculate the Euclidean distance
}

// Function to sort vertices based on proximity (nearest neighbor algorithm)
Eigen::MatrixXd sortVerticesByProximity(const Eigen::MatrixXd& V) {
    int numVertices = V.rows();
    Eigen::MatrixXd sortedVertices(numVertices, V.cols()); // To store the sorted vertices
    std::vector<bool> visited(numVertices, false);

    // Start with the first vertex
    int currentIndex = 0;
    sortedVertices.row(0) = V.row(currentIndex); // Add the first vertex to the sorted list
    visited[currentIndex] = true;

    // Continue selecting the nearest neighbor until all vertices are visited
    for (int i = 1; i < numVertices; ++i) {
        double minDistance = std::numeric_limits<double>::max();
        int nextIndex = -1;

        // Find the closest unvisited neighbor
        for (int j = 0; j < numVertices; ++j) {
            if (!visited[j]) {
                double distance = computeDistance(V.row(currentIndex), V.row(j)); // Full Euclidean distance
                if (distance < minDistance) {
                    minDistance = distance;
                    nextIndex = j;
                }
            }
        }

        // Move to the closest vertex if a valid next vertex is found
        if (nextIndex != -1) {
            currentIndex = nextIndex;
            sortedVertices.row(i) = V.row(currentIndex); // Add the next vertex to the sorted list
            visited[currentIndex] = true;
        }
        else {
            // This shouldn't happen if there are unvisited vertices
            std::cerr << "No unvisited neighbor found!" << std::endl;
            break;
        }
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
