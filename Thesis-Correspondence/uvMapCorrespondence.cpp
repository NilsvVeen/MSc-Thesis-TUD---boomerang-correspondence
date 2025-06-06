﻿#include <Eigen/Core>       // For Eigen::MatrixXd, Eigen::MatrixXi, and basic matrix operations
#include <Eigen/Dense>      // For additional dense matrix operations (e.g., decomposition)
#include <Eigen/Geometry>   // For geometric transformations if required (e.g., rotations)
#include <iostream>
#include <Eigen/Dense>
#include <vector>

#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>
#include <polyscope/surface_mesh.h>
#include "stl_utils.h"

const double epsilon = 1e-8;



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
#include <polyscope/curve_network.h>
#include "file_utils.h"

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

    std::cout << "Starting findConnectedBorder with B1 size: " << B1.rows() << " x " << B1.cols() << "\n";

    // Iterate over all pairs of consecutive boundary points (with wraparound for last case)
    for (int i = 0; i < B1.rows(); ++i) {
        int idxA = i;  // Boundary vertex a
        int idxB = (i + 1) % B1.rows();  // Boundary vertex b, with wraparound

        // Find corresponding indices of boundary vertices in V1
        Eigen::RowVector3d pointA = B1.row(idxA);
        Eigen::RowVector3d pointB = B1.row(idxB);

        //std::cout << "Processing boundary segment " << idxA << " -> " << idxB << "\n";
        //std::cout << "Boundary points: A = " << pointA << ", B = " << pointB << "\n";

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

        //std::cout << "Nearest vertices in V1: vIdxA = " << vIdxA << ", vIdxB = " << vIdxB << "\n";
        //std::cout << "Nearest vertices in V1: vA = " << V1.row(vIdxA) << ", vIdxB = " << V1.row(vIdxB) << "\n";
        //std::cout << "Min distances: minDistA = " << minDistA << ", minDistB = " << minDistB << "\n";

        if (vIdxA == -1 || vIdxB == -1) {
            std::cerr << "Error: Could not find nearest vertices for boundary points!\n";
            continue;
        }

        // Find the shortest path from vIdxA to vIdxB in V1
        std::vector<int> path = dijkstra(V1, F1, vIdxA, vIdxB);
        //std::cout << "Path from vIdxA to vIdxB: ";
        //for (int idx : path) {
        //    std::cout << idx << " ";
        //}
        //std::cout << "\n";

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

    std::cout << "Completed findConnectedBorder. Total unique vertices: " << uniqueVerticesOrdered.size() << "\n";
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

            int borderVertexCount = 0;
            for (int j = 0; j < 3; ++j) { // Loop through all vertices of the face
                Eigen::RowVector3d vertex = V1.row(F1(i, j));
                for (int k = 0; k < connectedBorder.rows(); ++k) {
                    if (vertex.isApprox(connectedBorder.row(k))) {
                        ++borderVertexCount;
                        if (borderVertexCount >= 2) {
                            break;  // If we already found 2 vertices, no need to check further
                        }
                    }
                }
                if (borderVertexCount >= 2) break;
            }
            if (borderVertexCount >= 2) {
                visited[i] = true;
                sideA.push_back(i); // Include in sideA since it touches the border
                continue;
            }


            //// If the neighboring face touches the border, we mark it as visited but don't process further
            //if (faceTouchesBorder(i)) {
            //    visited[i] = true;
            //    sideA.push_back(i); // Include in sideA since it touches the border
            //    continue;
            //}

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

void visualizeResults(
    const std::vector<Eigen::RowVector3d>& pointCloudA,
    const std::vector<Eigen::RowVector3d>& pointCloudA_skipped,
    const std::vector<Eigen::RowVector3d>& pointCloudC,
    const Eigen::MatrixXi& F1_new) {

    polyscope::init();

    // SECTION 1: Register Point Clouds
    {
        auto* pc1 = polyscope::registerPointCloud("P1 Original", pointCloudA);
        auto* pc1_skipped = polyscope::registerPointCloud("P1 Skipped", pointCloudA_skipped);
        auto* pc3_interpolated = polyscope::registerPointCloud("P3 Interpolated", pointCloudC);

        // Toggle visibility for point clouds
        pc1->setEnabled(true);                // Enable P1 point cloud by default
        pc1_skipped->setEnabled(false);       // Disable skipped points by default
        pc3_interpolated->setEnabled(true);   // Enable P3 interpolated points by default
    }

    // SECTION 2: Register Surface Meshes
    {
        polyscope::registerSurfaceMesh("M1 New (Original Mesh)", pointCloudA, F1_new);
        polyscope::registerSurfaceMesh("M2 New (Interpolated Mesh)", pointCloudC, F1_new);
    }

    // SECTION 3: Concatenate Point Clouds A and C
    std::vector<Eigen::RowVector3d> points = pointCloudA;
    points.insert(points.end(), pointCloudC.begin(), pointCloudC.end());

    // SECTION 4: Create Edges for Corresponding Points
    std::vector<std::array<int, 2>> edges;
    int numPoints = pointCloudA.size();  // Assume pointCloudA and pointCloudC have the same size
    for (int i = 0; i < numPoints; ++i) {
        edges.push_back({ i, numPoints + i });  // Corresponding points: A[i] -> C[i]
    }

    // SECTION 5: Split Curve Network into Smaller Parts
    int numParts = 100;  // You can adjust the number of parts based on your needs
    int partSize = numPoints / numParts;

    for (int part = 0; part < numParts; ++part) {
        // Determine the range of indices for this part
        int startIdx = part * partSize;
        int endIdx = (part == numParts - 1) ? numPoints : (part + 1) * partSize;

        // Extract the edges for this part
        std::vector<std::array<int, 2>> partEdges(edges.begin() + startIdx, edges.begin() + endIdx);

        // Register the curve network for this part
        auto* curveNetwork = polyscope::registerCurveNetwork("Corresponding Lines Part " + std::to_string(part + 1), points, partEdges);

        curveNetwork->setEnabled(false);

        // SECTION 6: Apply Rainbow Color Mapping for the Part
        Eigen::VectorXd colorValues(endIdx - startIdx);  // One value per edge for color mapping
        for (int i = startIdx; i < endIdx; ++i) {
            colorValues[i - startIdx] = static_cast<double>(i) / (numPoints - 1);  // Normalize the indices to [0, 1]
        }

        curveNetwork->addEdgeColorQuantity("Rainbow", colorValues);

        // Set the radius of the curve network
        curveNetwork->setRadius(0.005);
    }

    // Show Polyscope GUI
    polyscope::show();
}







double pointToTriangleDistance(
    const Eigen::RowVector2d& p,
    const Eigen::RowVector2d& A,
    const Eigen::RowVector2d& B,
    const Eigen::RowVector2d& C) {
    // Compute distances to the edges and vertices
    double d1 = (p - A).norm();
    double d2 = (p - B).norm();
    double d3 = (p - C).norm();
    return std::min({ d1, d2, d3 });
}

Eigen::RowVector2d projectPointOntoTriangle(
    const Eigen::RowVector2d& p,
    const Eigen::RowVector2d& A,
    const Eigen::RowVector2d& B,
    const Eigen::RowVector2d& C) {
    // For simplicity, return the closest vertex (full projection is more complex)
    double d1 = (p - A).squaredNorm();
    double d2 = (p - B).squaredNorm();
    double d3 = (p - C).squaredNorm();

    if (d1 < d2 && d1 < d3) return A;
    if (d2 < d1 && d2 < d3) return B;
    return C;
}




void UVToCorrespondence(
    const Eigen::MatrixXd& V1,  // Mesh 1 vertices
    const Eigen::MatrixXi& F1,  // Mesh 1 faces
    const Eigen::MatrixXd& B1,  // Boundary 1 vertices
    const Eigen::MatrixXd& UV1, // UV map of mesh 1

    const Eigen::MatrixXd& V2,  // Mesh 2 vertices
    const Eigen::MatrixXi& F2,  // Mesh 2 faces
    const Eigen::MatrixXd& B2,  // Boundary 2 vertices
    const Eigen::MatrixXd& UV2, // UV map of mesh 2

    const std::string correspondence3dMatched
) {
    // --- Settings and I/O
    Eigen::MatrixXd connectedBorder, connectedBorder2;
    bool readExisting = false;

    if (readExisting) {
        connectedBorder = readVerticesFromPLY(correspondence3dMatched + "/borders.obj");
        connectedBorder2 = readVerticesFromPLY(correspondence3dMatched + "/borders2.obj");
    }
    else {
        connectedBorder = findConnectedBorder(V1, F1, B1);
        writeVerticesToPLY(correspondence3dMatched + "/borders.obj", connectedBorder);

        connectedBorder2 = findConnectedBorder(V2, F2, B2);
        writeVerticesToPLY(correspondence3dMatched + "/borders2.obj", connectedBorder2);
    }

    // --- Side classification
    auto [sideA, sideB] = classifyFacesByBorder(V1, F1, connectedBorder);
    auto [sideA2, sideB2] = classifyFacesByBorder(V2, F2, connectedBorder2);

    std::cout << "Mesh 1 - Side A: " << sideA.size() << ", Side B: " << sideB.size() << std::endl;
    std::cout << "Mesh 2 - Side A: " << sideA2.size() << ", Side B: " << sideB2.size() << std::endl;

    // --- Store classifications globally/static for GUI toggle
    static std::vector<int> sideA_stored = sideA;
    static std::vector<int> sideB_stored = sideB;
    static std::vector<int> sideA2_stored = sideA2;
    static std::vector<int> sideB2_stored = sideB2;
    static bool flipped = false;

    // --- Polyscope setup
    polyscope::init();

    // Register border and boundary points
    polyscope::registerPointCloud("Connected Border 1", connectedBorder);
    polyscope::registerPointCloud("Original Boundary 1", B1);
    polyscope::registerPointCloud("Connected Border 2", connectedBorder2);
    polyscope::registerPointCloud("Original Boundary 2", B2);

    // Register surface meshes
    polyscope::registerSurfaceMesh("M111", V1, F1);
    polyscope::registerSurfaceMesh("M222", V2, F2);

    // --- Scalars
    static Eigen::VectorXd faceColors(F1.rows());
    static Eigen::VectorXd faceColors2(F2.rows());

    auto updateFaceClassificationDisplay = [&]() {
        faceColors.setConstant(-1);
        faceColors2.setConstant(-1);

        const auto& currentA1 = sideA_stored;
        const auto& currentB1 = sideB_stored;
        const auto& currentA2 = flipped ? sideB2_stored : sideA2_stored;
        const auto& currentB2 = flipped ? sideA2_stored : sideB2_stored;

        for (int faceIdx : currentA1) faceColors(faceIdx) = 0;
        for (int faceIdx : currentB1) faceColors(faceIdx) = 1;
        for (int faceIdx : currentA2) faceColors2(faceIdx) = 0;
        for (int faceIdx : currentB2) faceColors2(faceIdx) = 1;

        polyscope::getSurfaceMesh("M111")->addFaceScalarQuantity("Side Classification", faceColors, polyscope::DataType::SYMMETRIC)->setEnabled(true);
        polyscope::getSurfaceMesh("M222")->addFaceScalarQuantity("Side Classification2", faceColors2, polyscope::DataType::SYMMETRIC)->setEnabled(true);
    };

    // First display
    updateFaceClassificationDisplay();

    // Add a GUI button to flip
    polyscope::state::userCallback = [&]() {
        if (ImGui::Button("Swap Sides")) {
            flipped = !flipped;
            updateFaceClassificationDisplay();
        }
    };

    // Launch GUI
    polyscope::show();

    if (flipped) {
        sideA2 = sideB2_stored;
        sideB2 = sideA2_stored;
    }



    std::vector<Eigen::RowVector3d> pointCloudA;          // Holds 3D points from V1
    std::vector<Eigen::RowVector3d> pointCloudA_skipped;  // Holds skipped points from V1
    std::vector<Eigen::RowVector3d> pointCloudC;          // Holds interpolated points in V2
    std::vector<Eigen::RowVector2d> pointCloudC_UV;          // Holds interpolated points in V2



    std::vector<int> indicesSkipped;
    for (int i = 0; i < V1.rows(); ++i) {
        Eigen::RowVector3d a = V1.row(i);   // Vertex in 3D
        Eigen::RowVector2d b = UV1.row(i); // Corresponding UV coordinate

        // Check which side `a` belongs to
        bool isInSideA = false, isInSideB = false;
        for (int f : sideA) {
            if ((F1.row(f).array() == i).any()) {  // If vertex `i` is part of a face in `sideA`
                isInSideA = true;
                break;
            }
        }
        if (!isInSideA) {
            for (int f : sideB) {
                if ((F1.row(f).array() == i).any()) {  // If vertex `i` is part of a face in `sideB`
                    isInSideB = true;
                    break;
                }
            }
        }

        if (!isInSideA && !isInSideB) {
            std::cerr << "Error: Vertex " << i << " does not belong to any side.\n";
            pointCloudA_skipped.push_back(a);
            continue;
        }

        bool added = false;

        // matching side cases
        // Find the face in UV2 that contains `b`
        for (int j = 0; j < F2.rows(); ++j) {
            // Get the UV coordinates of the face vertices in UV2
            Eigen::RowVector2d A = UV2.row(F2(j, 0));
            Eigen::RowVector2d B = UV2.row(F2(j, 1));
            Eigen::RowVector2d C = UV2.row(F2(j, 2));

            // Check if `b` lies in the UV triangle
            if (isPointInTriangle(b, A, B, C)) {
                // Check if this face belongs to the correct side
                bool isInSideA2 = std::find(sideA2.begin(), sideA2.end(), j) != sideA2.end();
                bool isInSideB2 = std::find(sideB2.begin(), sideB2.end(), j) != sideB2.end();

                if ((isInSideA && !isInSideA2) || (isInSideB && !isInSideB2)) {
                    // Skip the triangle if its on the wrong side
                    continue;
                }

                // Compute barycentric coordinates
                Eigen::Matrix2d T;
                T.col(0) = B - A;
                T.col(1) = C - A;
                Eigen::RowVector2d v = b - A;

                Eigen::Vector2d barycentric2D;
                try {
                    barycentric2D = T.inverse() * v.transpose();
                }
                catch (...) {

                    std::cout << "???????WARNING" << std::endl;
                    continue; // Skip this triangle if there's an error
                }

                double lambda1 = 1.0 - barycentric2D.sum();
                double lambda2 = barycentric2D[0];
                double lambda3 = barycentric2D[1];

                // Validate barycentric coordinates
                if (lambda1 < -epsilon || lambda2 < -epsilon || lambda3 < -epsilon ||
                    lambda1 > 1.0 + epsilon || lambda2 > 1.0 + epsilon || lambda3 > 1.0 + epsilon) {
                    std::cout << "Invalid Barycentric Coordinates for Vertex " << i << ":\n";
                    std::cout << "  lambda1 = " << lambda1 << ", lambda2 = " << lambda2 << ", lambda3 = " << lambda3 << "\n";

                    continue; // Skip invalid coordinates
                }

                // Interpolate 3D position in Mesh 2
                Eigen::RowVector3d interpolatedPoint =
                    lambda1 * V2.row(F2(j, 0)) +
                    lambda2 * V2.row(F2(j, 1)) +
                    lambda3 * V2.row(F2(j, 2));

                // Additional sanity check for the interpolated point
                if (interpolatedPoint.hasNaN() || interpolatedPoint.norm() > 1e6) {

                    std::cout << "???????out???" << std::endl;

                    continue; // Skip out-of-range points
                }

                // Add the vertex and its interpolated point
                pointCloudA.push_back(a);
                pointCloudC.push_back(interpolatedPoint);
                added = true;
                break; // Move to the next vertex in V1
            }
        }









        // not matching side, but needed for near boundary cases
        if (!added) {
            // Find the face in UV2 that contains `b`
            for (int j = 0; j < F2.rows(); ++j) {
                // Get the UV coordinates of the face vertices in UV2
                Eigen::RowVector2d A = UV2.row(F2(j, 0));
                Eigen::RowVector2d B = UV2.row(F2(j, 1));
                Eigen::RowVector2d C = UV2.row(F2(j, 2));

                // Check if `b` lies in the UV triangle
                if (isPointInTriangle(b, A, B, C)) {

                    // Compute barycentric coordinates
                    Eigen::Matrix2d T;
                    T.col(0) = B - A;
                    T.col(1) = C - A;
                    Eigen::RowVector2d v = b - A;

                    Eigen::Vector2d barycentric2D;
                    try {
                        barycentric2D = T.inverse() * v.transpose();
                    }
                    catch (...) {

                        std::cout << "???????WARNING" << std::endl;
                        continue; // Skip this triangle if there's an error
                    }

                    double lambda1 = 1.0 - barycentric2D.sum();
                    double lambda2 = barycentric2D[0];
                    double lambda3 = barycentric2D[1];

                    // Validate barycentric coordinates
                    if (lambda1 < -epsilon || lambda2 < -epsilon || lambda3 < -epsilon ||
                        lambda1 > 1.0 + epsilon || lambda2 > 1.0 + epsilon || lambda3 > 1.0 + epsilon) {
                        std::cout << "Invalid Barycentric Coordinates for Vertex " << i << ":\n";
                        std::cout << "  lambda1 = " << lambda1 << ", lambda2 = " << lambda2 << ", lambda3 = " << lambda3 << "\n";

                        continue; // Skip invalid coordinates
                    }

                    // Interpolate 3D position in Mesh 2
                    Eigen::RowVector3d interpolatedPoint =
                        lambda1 * V2.row(F2(j, 0)) +
                        lambda2 * V2.row(F2(j, 1)) +
                        lambda3 * V2.row(F2(j, 2));

                    // Additional sanity check for the interpolated point
                    if (interpolatedPoint.hasNaN() || interpolatedPoint.norm() > 1e6) {

                        std::cout << "???????out???" << std::endl;

                        continue; // Skip out-of-range points
                    }

                    // Add the vertex and its interpolated point
                    pointCloudA.push_back(a);
                    pointCloudC.push_back(interpolatedPoint);
                    added = true;
                    break; // Move to the next vertex in V1
                }
            }
        }

        // todo
        // in this case it is no triangles, but outside the uv map. 
        if (!added) {
            // Find the nearest face in UV2 to `b`
            double minDistance = std::numeric_limits<double>::max();
            int nearestFaceIdx = -1;

            for (int j = 0; j < F2.rows(); ++j) {
                // Get the UV coordinates of the face vertices in UV2
                Eigen::RowVector2d A = UV2.row(F2(j, 0));
                Eigen::RowVector2d B = UV2.row(F2(j, 1));
                Eigen::RowVector2d C = UV2.row(F2(j, 2));

                // Compute the distance from `b` to the face
                double distance = pointToTriangleDistance(b, A, B, C);

                if (distance < minDistance) {
                    minDistance = distance;
                    nearestFaceIdx = j;
                }
            }

            if (nearestFaceIdx != -1) {
                // Project `b` onto the nearest face
                Eigen::RowVector2d A = UV2.row(F2(nearestFaceIdx, 0));
                Eigen::RowVector2d B = UV2.row(F2(nearestFaceIdx, 1));
                Eigen::RowVector2d C = UV2.row(F2(nearestFaceIdx, 2));

                Eigen::RowVector2d projectedPoint = projectPointOntoTriangle(b, A, B, C);

                // Compute barycentric coordinates for the projected point
                Eigen::Matrix2d T;
                T.col(0) = B - A;
                T.col(1) = C - A;
                Eigen::RowVector2d v = projectedPoint - A;

                Eigen::Vector2d barycentric2D;
                try {
                    barycentric2D = T.inverse() * v.transpose();
                }
                catch (...) {
                    std::cerr << "WARNING: Barycentric coordinate calculation failed.\n";
                    return;
                }

                double lambda1 = 1.0 - barycentric2D.sum();
                double lambda2 = barycentric2D[0];
                double lambda3 = barycentric2D[1];

                // Interpolate the 3D position in V2
                Eigen::RowVector3d interpolatedPoint =
                    lambda1 * V2.row(F2(nearestFaceIdx, 0)) +
                    lambda2 * V2.row(F2(nearestFaceIdx, 1)) +
                    lambda3 * V2.row(F2(nearestFaceIdx, 2));

                // Add the vertex and its interpolated point
                added = true;
                pointCloudA.push_back(a);
                pointCloudC.push_back(interpolatedPoint);
            }
            else {
                // If no face is found (extremely rare), log and skip the point
                std::cerr << "WARNING: Could not find a nearest face for vertex " << i << ".\n";
            }
        }


        if (!added) {
            pointCloudC_UV.push_back(b);
            pointCloudA_skipped.push_back(a);
            indicesSkipped.push_back(i);
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

    // Print summary of skipped points
    std::cout << "\nSkipped vertices summary:\n";
    std::cout << "  original skipped: " << pointCloudA_skipped.size() << "\n";
    std::cout << "  A: " << pointCloudA.size() << "\n";
    std::cout << "  C: " << pointCloudC.size() << "\n";


 polyscope::init();
    // Visualize in Polyscope
polyscope::registerSurfaceMesh("Mesh", V1, F1); // Register the first mesh
polyscope::registerSurfaceMesh("Mesh 2", V2, F2); // Register the second mesh
polyscope::registerSurfaceMesh2D("UV M`1", UV1, F1); // Register the second mesh
polyscope::registerSurfaceMesh2D("UV M`2", UV2, F2); // Register the second mesh
polyscope::registerPointCloud2D("skipped points", pointCloudC_UV); // Register the second mesh


// Show Polyscope (blocks execution until window is closed)
polyscope::show();
polyscope::removeAllGroups();
polyscope::removeAllStructures();





    std::cout << " length: " << F1.rows() << std::endl;
    auto VF = filterMeshUsingPoints(V1, F1, pointCloudA_skipped);
    auto F1_new = VF.second;
    auto V1_new = VF.first;
    std::cout << " length2 : " << F1_new.rows() << std::endl;

    visualizeResults(
        pointCloudA,
        pointCloudA_skipped,
        pointCloudC,
        F1_new);


    polyscope::show();

    // Debug output
    std::cout << "Registered point clouds:" << std::endl;
    std::cout << "PointCloudA size: " << pointCloudA.size() << std::endl;
    std::cout << "SKIPPED: " << pointCloudA_skipped.size() << std::endl;
    std::cout << "PointCloudC size: " << pointCloudC.size() << std::endl;

    //createDirectory(correspondence3dMatched);
    //clearDirectory(correspondence3dMatched);


    Eigen::MatrixXd pointCloudMatrixA(pointCloudA.size(), 3);  // Create a matrix with the same number of rows as the vector
    for (size_t i = 0; i < pointCloudA.size(); ++i) {
        pointCloudMatrixA.row(i) = pointCloudA[i];  // Copy each RowVector3d into a row of the matrix
    }

    Eigen::MatrixXd pointCloudMatrixC(pointCloudC.size(), 3);  // Create a matrix with the same number of rows as the vector
    for (size_t i = 0; i < pointCloudC.size(); ++i) {
        pointCloudMatrixC.row(i) = pointCloudC[i];  // Copy each RowVector3d into a row of the matrix
    }


    saveMeshToFile(correspondence3dMatched + "/M1.obj", pointCloudMatrixA, F1_new);
    saveMeshToFile(correspondence3dMatched + "/M2.obj", pointCloudMatrixC, F1_new);
   

}

