#include <Eigen/Core>       // For Eigen::MatrixXd, Eigen::MatrixXi, and basic matrix operations
#include <Eigen/Dense>      // For additional dense matrix operations (e.g., decomposition)
#include <Eigen/Geometry>   // For geometric transformations if required (e.g., rotations)
#include <iostream>
#include <Eigen/Dense>
#include <vector>

#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>
#include <polyscope/surface_mesh.h>
#include "stl_utils.h"





// Function to check if a point is inside a triangle in UV space
bool isPointInTriangle(const Eigen::RowVector2d& point, const Eigen::RowVector2d& A,
    const Eigen::RowVector2d& B, const Eigen::RowVector2d& C) {
    Eigen::RowVector2d v0 = C - A;
    Eigen::RowVector2d v1 = B - A;
    Eigen::RowVector2d v2 = point - A;

    double dot00 = v0.dot(v0);
    double dot01 = v0.dot(v1);
    double dot02 = v0.dot(v2);
    double dot11 = v1.dot(v1);
    double dot12 = v1.dot(v2);

    double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    return (u >= 0) && (v >= 0) && (u + v <= 1);
}


std::pair<Eigen::MatrixXd, Eigen::MatrixXi> filterMeshUsingPoints(
    const Eigen::MatrixXd& V1,                  // Original vertices
    const Eigen::MatrixXi& F1,                 // Original faces
    const std::vector<Eigen::RowVector3d>& pointstoremove // Points to be removed
) {
    // Step 1: Convert faces to point-based format
    std::vector<std::array<Eigen::RowVector3d, 3>> pointBasedFaces;
    for (int i = 0; i < F1.rows(); ++i) {
        pointBasedFaces.push_back({
            V1.row(F1(i, 0)),
            V1.row(F1(i, 1)),
            V1.row(F1(i, 2))
            });
    }

    // Hash function for Eigen::RowVector3d
    struct hashRowVector3d {
        std::size_t operator()(const Eigen::RowVector3d& vec) const {
            std::size_t h1 = std::hash<double>()(vec.x());
            std::size_t h2 = std::hash<double>()(vec.y());
            std::size_t h3 = std::hash<double>()(vec.z());
            return h1 ^ (h2 << 1) ^ (h3 << 2); // Combine hashes
        }
    };

    // Equality function for Eigen::RowVector3d
    struct equalRowVector3d {
        bool operator()(const Eigen::RowVector3d& lhs, const Eigen::RowVector3d& rhs) const {
            return lhs.isApprox(rhs); // Approximate equality comparison
        }
    };

    // Step 2: Remove points from the vertex list
    std::unordered_set<int> removedIndices;
    Eigen::MatrixXd newV(V1.rows(), 3);
    int newIndex = 0;
    std::unordered_map<Eigen::RowVector3d, int, hashRowVector3d> pointToIndex; // Map for fast lookup
    for (int i = 0; i < V1.rows(); ++i) {
        bool toRemove = false;
        for (const auto& point : pointstoremove) {
            if (point.isApprox(V1.row(i))) {
                toRemove = true;
                removedIndices.insert(i);
                break;
            }
        }

        if (!toRemove) {
            newV.row(newIndex) = V1.row(i);
            pointToIndex[V1.row(i)] = newIndex;
            ++newIndex;
        }
    }
    newV.conservativeResize(newIndex, Eigen::NoChange);

    // Step 3: Rebuild faces with updated indices
    std::vector<Eigen::RowVector3i> newFaces;
    for (const auto& face : pointBasedFaces) {
        bool isValidFace = true;
        Eigen::RowVector3i newFace;

        for (int j = 0; j < 3; ++j) {
            const auto& vertex = face[j];
            if (pointToIndex.find(vertex) == pointToIndex.end()) {
                isValidFace = false; // One of the vertices was removed
                break;
            }
            newFace[j] = pointToIndex[vertex];
        }

        if (isValidFace) {
            newFaces.push_back(newFace);
        }
    }

    // Convert the face list to an Eigen matrix
    Eigen::MatrixXi newF(newFaces.size(), 3);
    for (size_t i = 0; i < newFaces.size(); ++i) {
        newF.row(i) = newFaces[i];
    }

    return { newV, newF };
}



#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <limits>

using namespace std;

// Euclidean distance between two points
double euclideanDistance(const Eigen::RowVector3d& p1, const Eigen::RowVector3d& p2) {
    return (p1 - p2).norm();
}

// Build adjacency list based on the mesh faces
unordered_map<int, unordered_set<int>> buildAdjacencyList(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    unordered_map<int, unordered_set<int>> adjList;
    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = j + 1; k < 3; ++k) {
                adjList[F(i, j)].insert(F(i, k));
                adjList[F(i, k)].insert(F(i, j));
            }
        }
    }
    return adjList;
}

// Dijkstra's algorithm to find the shortest path between two vertices
vector<int> dijkstra(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int startIdx, int endIdx) {
    int numVertices = V.rows();
    vector<double> dist(numVertices, std::numeric_limits<double>::infinity());
    vector<int> prev(numVertices, -1);
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;

    dist[startIdx] = 0.0;
    pq.push({ 0.0, startIdx });

    // Build adjacency list from faces
    unordered_map<int, unordered_set<int>> adjList = buildAdjacencyList(V, F);

    // Dijkstra's algorithm
    unordered_set<int> visited; // To avoid revisiting vertices
    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();

        if (visited.find(u) != visited.end()) {
            continue; // Skip if already visited
        }

        visited.insert(u); // Mark as visited

        if (u == endIdx) break; // Short-circuit if we reached the target

        for (int v : adjList[u]) {
            if (visited.find(v) != visited.end()) {
                continue; // Skip if already visited
            }
            double weight = euclideanDistance(V.row(u), V.row(v));
            if (dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                prev[v] = u;
                pq.push({ dist[v], v });
            }
        }
    }

    // Reconstruct the path from start to end
    vector<int> path;
    for (int u = endIdx; u != -1; u = prev[u]) {
        path.push_back(u);
    }

    reverse(path.begin(), path.end());
    return path;
}



Eigen::MatrixXd findConnectedBorder(const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F1, const Eigen::MatrixXd& B1) {
    // List to store unique vertices in order
    std::vector<Eigen::RowVector3d> uniqueVerticesOrdered;

    // Function to add unique vertices while preserving order
    auto addUniqueVertex = [&uniqueVerticesOrdered](const Eigen::RowVector3d& vertex) {
        // Check if the vertex already exists in the list
        bool isDuplicate = false;
        for (const auto& v : uniqueVerticesOrdered) {
            if (v.isApprox(vertex)) {  // Use isApprox for floating-point comparison
                isDuplicate = true;
                break;
            }
        }
        // If not duplicate, add it to the list
        if (!isDuplicate) {
            uniqueVerticesOrdered.push_back(vertex);
        }
    };

    // Iterate over all pairs of consecutive boundary points (with wraparound for last case)
    for (int i = 0; i < B1.rows(); ++i) {
        int idxA = i;  // Boundary vertex a
        int idxB = (i + 1) % B1.rows();  // Boundary vertex b, with wraparound

        // Find corresponding indices of boundary vertices in V1
        Eigen::RowVector3d pointA = B1.row(idxA);
        Eigen::RowVector3d pointB = B1.row(idxB);

        // Match B1.row(idxA) and B1.row(idxB) to the nearest vertices in V1
        int vIdxA = -1, vIdxB = -1;
        double minDistA = std::numeric_limits<double>::infinity(), minDistB = std::numeric_limits<double>::infinity();

        for (int j = 0; j < V1.rows(); ++j) {
            double distA = (pointA - V1.row(j)).norm();
            double distB = (pointB - V1.row(j)).norm();

            if (distA < minDistA) {
                minDistA = distA;
                vIdxA = j;
            }
            if (distB < minDistB) {
                minDistB = distB;
                vIdxB = j;
            }
        }

        // Find the shortest path from vIdxA to vIdxB in V1
        std::vector<int> path = dijkstra(V1, F1, vIdxA, vIdxB);

        // Add all vertices in the path to the list, ensuring uniqueness
        for (int idx : path) {
            addUniqueVertex(V1.row(idx));
        }
    }

    // Convert uniqueVerticesOrdered to a matrix
    Eigen::MatrixXd uniqueVerticesMat(uniqueVerticesOrdered.size(), 3);
    for (size_t i = 0; i < uniqueVerticesOrdered.size(); ++i) {
        uniqueVerticesMat.row(i) = uniqueVerticesOrdered[i];
    }

    return uniqueVerticesMat;
}





std::pair<std::vector<int>, std::vector<int>> classifyFacesByBorder(
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& F1,
    const Eigen::MatrixXd& connectedBorder)
{
    // Data structures to store classified faces
    std::vector<int> sideA, sideB;
    std::vector<bool> visited(F1.rows(), false);

    // Helper function to determine if a face contains a vertex from the connected border
    auto faceTouchesBorder = [&](int faceIdx) {
        for (int i = 0; i < 3; ++i) {
            Eigen::RowVector3d vertex = V1.row(F1(faceIdx, i));
            for (int j = 0; j < connectedBorder.rows(); ++j) {
                if (vertex.isApprox(connectedBorder.row(j))) {
                    return true;
                }
            }
        }
        return false;
    };

    // Helper function to check if two faces share an edge
    auto facesShareEdge = [&](int f1, int f2) {
        int sharedVertices = 0;
        for (int vi = 0; vi < 3; ++vi) {
            for (int vj = 0; vj < 3; ++vj) {
                if (F1(f1, vi) == F1(f2, vj)) {
                    ++sharedVertices;
                }
            }
        }
        return sharedVertices >= 2; // Shared edge means at least 2 vertices in common
    };

    // Start classification from a seed face (on side A) that touches the border
    int seedFace = -1;
    for (int i = 0; i < F1.rows(); ++i) {
        if (faceTouchesBorder(i)) {
            seedFace = i;
            break;
        }
    }
    if (seedFace == -1) {
        throw std::runtime_error("No seed face found adjacent to the border!");
    }

    // Flood-fill to classify the faces
    std::queue<int> faceQueue;
    faceQueue.push(seedFace);
    visited[seedFace] = true;
    sideA.push_back(seedFace);

    while (!faceQueue.empty()) {
        int currentFace = faceQueue.front();
        faceQueue.pop();

        // Check all neighbors of the current face
        for (int i = 0; i < F1.rows(); ++i) {
            if (visited[i] || !facesShareEdge(currentFace, i)) continue; // Skip visited or unconnected faces

            // If the neighboring face touches the border, we mark it as visited but don't process further
            if (faceTouchesBorder(i)) {
                visited[i] = true;
                sideA.push_back(i); // Include in sideA since it touches the border
                continue;
            }

            // Assign the face to the same side as the current face
            visited[i] = true;
            sideA.push_back(i);
            faceQueue.push(i);
        }
    }

    // Remaining unvisited faces are on the other side
    for (int i = 0; i < F1.rows(); ++i) {
        if (!visited[i]) {
            sideB.push_back(i);
        }
    }

    return { sideA, sideB };
}




void UVToCorrespondence(
    const Eigen::MatrixXd& V1,  // Mesh 1 vertices       n x 3
    const Eigen::MatrixXi& F1, // Mesh 1 faces          m x 3
    const Eigen::MatrixXd& B1, // Boundary 1 vertices   p x 3
    const Eigen::MatrixXd& UV1, // UV map of mesh 1     n x 2

    const Eigen::MatrixXd& V2,  // Mesh 2 vertices       s x 3
    const Eigen::MatrixXi& F2,  // Mesh 2 faces          q x 3
    const Eigen::MatrixXd& B2,  // Boundary 2 vertices   p x 3
    const Eigen::MatrixXd& UV2  // UV map of mesh 2     s x 2
) {

//Eigen::MatrixXd connectedBorder = findConnectedBorder(V1, F1, B1);
//
//writeVerticesToPLY("borders.obj", connectedBorder);

Eigen::MatrixXd connectedBorder = readVerticesFromPLY("borders.obj");

auto [sideA, sideB] = classifyFacesByBorder(V1, F1, connectedBorder);

// Print the face indices for each side
std::cout << "Faces on Side A: " << sideA.size() << std::endl;
std::cout << "Faces on Side B: " << sideB.size() << std::endl;

// Create a scalar field to color faces based on their classification
Eigen::VectorXd faceColors(F1.rows());
faceColors.setConstant(-1); // Default value for unclassified faces (optional)

// Assign colors to faces in sideA and sideB
for (int faceIdx : sideA) {
    faceColors(faceIdx) = 0; // Color for Side A
}
for (int faceIdx : sideB) {
    faceColors(faceIdx) = 1; // Color for Side B
}

polyscope::init();

// Register connected border as a point cloud
polyscope::registerPointCloud("Unique Vertices", connectedBorder);
// Register the original boundary vertices
polyscope::registerPointCloud("Original Boundary Vertices", B1);

// Register the mesh with the face scalar quantity
polyscope::registerSurfaceMesh("M1", V1, F1)
    ->addFaceScalarQuantity("Side Classification", faceColors, polyscope::DataType::SYMMETRIC);

// Optionally register another mesh
polyscope::registerSurfaceMesh("M2", V2, F2);

polyscope::show();



    // Point clouds for registration
    std::vector<Eigen::RowVector3d> pointCloudA;          // Holds 3D points from V1
    std::vector<Eigen::RowVector3d> pointCloudA_skipped;  // Holds skipped points from V1
    std::vector<Eigen::RowVector3d> pointCloudC;          // Holds interpolated points in V2

    for (int i = 0; i < V1.rows(); ++i) {
        Eigen::RowVector3d a = V1.row(i);   // Vertex in 3D
        Eigen::RowVector2d b = UV1.row(i); // Corresponding UV coordinate

        bool added = false;
        // Find the face in UV2 that contains b
        for (int j = 0; j < F2.rows(); ++j) {
            // Get the UV coordinates of the face vertices in UV2
            Eigen::RowVector2d A = UV2.row(F2(j, 0));
            Eigen::RowVector2d B = UV2.row(F2(j, 1));
            Eigen::RowVector2d C = UV2.row(F2(j, 2));

            int count = 0;
            if (isPointInTriangle(b, A, B, C)) {
                Eigen::Matrix2d T;
                T.col(0) = B - A;
                T.col(1) = C - A;
                Eigen::RowVector2d v = b - A;

                Eigen::Vector2d barycentric2D;
                try {
                    barycentric2D = T.inverse() * v.transpose();
                }
                catch (...) {
                    continue; // Skip this triangle if there's an error
                }

                double lambda1 = 1.0 - barycentric2D.sum();
                double lambda2 = barycentric2D[0];
                double lambda3 = barycentric2D[1];

                // Validate barycentric coordinates
                if (lambda1 < 0.0 || lambda2 < 0.0 || lambda3 < 0.0 ||
                    lambda1 > 1.0 || lambda2 > 1.0 || lambda3 > 1.0) {
                    continue; // Skip invalid coordinates
                }

                // Interpolate 3D position in Mesh 2
                Eigen::RowVector3d interpolatedPoint =
                    lambda1 * V2.row(F2(j, 0)) +
                    lambda2 * V2.row(F2(j, 1)) +
                    lambda3 * V2.row(F2(j, 2));

                // Additional sanity check for the interpolated point
                if (interpolatedPoint.hasNaN() || interpolatedPoint.norm() > 1e6) {
                    continue; // Skip out-of-range points
                }
                count += 1;
                pointCloudA.push_back(a);
                pointCloudC.push_back(interpolatedPoint);
                added = true;
                break; // Move to the next vertex in V1
            }
        }
        if (!added) {
            pointCloudA_skipped.push_back(a);
        }

        // Progress bar logic
        int progressWidth = 50;
        double progress = static_cast<double>(i + 1) / V1.rows();
        int pos = static_cast<int>(progress * progressWidth);
        std::cout << "\r[";
        for (int k = 0; k < progressWidth; ++k) {
            if (k < pos) std::cout << "=";
            else if (k == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << std::fixed << std::setprecision(1)
            << (progress * 100.0) << "% completed" << std::flush;
    }

    std::cout << " length: " << F1.rows() << std::endl;
    auto VF = filterMeshUsingPoints(V1, F1, pointCloudA_skipped);
    auto F1_new = VF.second;
    auto V1_new = VF.first;
    std::cout << " length2 : " << F1_new.rows() << std::endl;


    // Register point clouds in Polyscope
    polyscope::registerPointCloud("p1", pointCloudA);
    polyscope::registerPointCloud("p1 skipped", pointCloudA_skipped);
    polyscope::registerPointCloud("p3 interpolated", pointCloudC);

    polyscope::registerSurfaceMesh("M1_new", pointCloudA, F1_new);
    polyscope::registerSurfaceMesh("M2_new", pointCloudC, F1_new);


    polyscope::show();

    // Debug output
    std::cout << "Registered point clouds:" << std::endl;
    std::cout << "PointCloudA size: " << pointCloudA.size() << std::endl;
    std::cout << "PointCloudC size: " << pointCloudC.size() << std::endl;
}

