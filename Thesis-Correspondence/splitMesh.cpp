#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <set>
#include <vector>
#include <unordered_map>
#include <stack>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/edge_lengths.h>
#include <igl/adjacency_matrix.h>
#include <iostream>
#include <Eigen/Dense>
#include <unordered_set>
#include <vector>
#include <iostream>

#include <Eigen/Dense>
#include <igl/connected_components.h>
#include <iostream>
#include <Eigen/Dense>
#include <igl/boundary_loop.h>
#include <unordered_set>
#include <vector>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <Eigen/Dense>


void removeVerticesWithTwoFacesAndBorderEdges(
    const Eigen::MatrixXd& MeshA_V,
    const Eigen::MatrixXi& MeshA_F,
    Eigen::MatrixXd& border_V, // Matrix of border vertices
    Eigen::MatrixXd& MeshB_V,
    Eigen::MatrixXi& MeshB_F,
    Eigen::MatrixXd& removed_V, // Matrix of removed vertices
    Eigen::MatrixXi& removed_F, // Matrix of removed faces
    Eigen::MatrixXd& border_V_new // Matrix of border vertices excluding removed ones
) {
    // Print input sizes
    std::cout << "Input MeshA_V size: " << MeshA_V.rows() << " vertices" << std::endl;
    std::cout << "Input MeshA_F size: " << MeshA_F.rows() << " faces" << std::endl;
    std::cout << "Input border_V size: " << border_V.rows() << " border vertices" << std::endl;

    // Set to store indices of border vertices in MeshA_V based on coordinate matching
    std::unordered_set<int> borderVertexIndices;
    double tolerance = 1e-4;

    // Match each vertex in border_V to a corresponding vertex in MeshA_V by coordinates
    for (int i = 0; i < border_V.rows(); ++i) {
        bool found = false;
        for (int j = 0; j < MeshA_V.rows(); ++j) {
            double distance = (MeshA_V.row(j) - border_V.row(i)).norm();
            if (distance < tolerance) {
                borderVertexIndices.insert(j);
                found = true;
                break;
            }
        }
        if (!found) {
            std::cout << "Warning: No close match found in MeshA_V for border vertex at "
                << border_V.row(i) << std::endl;
        }
    }

    // Set to keep track of vertices to remove based on the conditions
    std::unordered_set<int> verticesToRemove;
    std::vector<Eigen::Vector3i> removedFacesList;

    // Iterate over each border vertex index
    for (const int borderVertexIndex : borderVertexIndices) {
        std::vector<int> facesContainingVertex;

        // Collect all faces containing this border vertex
        for (int f = 0; f < MeshA_F.rows(); ++f) {
            for (int j = 0; j < 3; ++j) {
                if (MeshA_F(f, j) == borderVertexIndex) {
                    facesContainingVertex.push_back(f);
                    break;
                }
            }
        }

        // Apply the rules based on the number of faces containing this vertex
        int faceCount = facesContainingVertex.size();
        bool shouldRemove = false;

        // Rule based on unique border vertices count
        int uniqueBorderVertexCount = 0;
        std::unordered_set<int> uniqueBorderVertices;
        for (int face : facesContainingVertex) {
            for (int j = 0; j < 3; ++j) {
                int vertexIndex = MeshA_F(face, j);
                if (vertexIndex != borderVertexIndex && borderVertexIndices.count(vertexIndex) > 0) {
                    uniqueBorderVertices.insert(vertexIndex);
                }
            }
        }
        uniqueBorderVertexCount = uniqueBorderVertices.size();

        // Ensure that the sum of the unique vertices in faces containing border vertex x is exactly 4
        if (uniqueBorderVertexCount == 4) {
            // Additional conditions based on the face count
            if (faceCount == 2) {
                // Rule 1: All vertices in both faces must be border vertices
                shouldRemove = true;
                for (int face : facesContainingVertex) {
                    for (int j = 0; j < 3; ++j) {
                        int vertexIndex = MeshA_F(face, j);
                        if (vertexIndex != borderVertexIndex && borderVertexIndices.find(vertexIndex) == borderVertexIndices.end()) {
                            shouldRemove = false;
                            break;
                        }
                    }
                    if (!shouldRemove) break;
                }
            }
            else if (faceCount == 3) {
                // Rule 2: One face should have 3 border vertices, two should have 2
                int facesWithThreeBorderVertices = 0;
                int facesWithTwoBorderVertices = 0;

                for (int face : facesContainingVertex) {
                    int borderVertexCount = 0;
                    for (int j = 0; j < 3; ++j) {
                        int vertexIndex = MeshA_F(face, j);
                        if (vertexIndex == borderVertexIndex || borderVertexIndices.count(vertexIndex) > 0) {
                            borderVertexCount++;
                        }
                    }
                    if (borderVertexCount == 3) facesWithThreeBorderVertices++;
                    if (borderVertexCount == 2) facesWithTwoBorderVertices++;
                }
                shouldRemove = (facesWithThreeBorderVertices == 1 && facesWithTwoBorderVertices == 2);
            }
            else if (faceCount == 4) {
                // todo fix this add conditions
                shouldRemove = true;
            }
        }

        // If all conditions are met, mark this vertex for removal and add its faces to removed_F
        if (shouldRemove) {
            verticesToRemove.insert(borderVertexIndex);
            for (int face : facesContainingVertex) {
                removedFacesList.emplace_back(MeshA_F.row(face));
            }
        }
    }

    // Update removed_F with the faces associated with removed vertices
    removed_F.resize(removedFacesList.size(), 3);
    for (size_t i = 0; i < removedFacesList.size(); ++i) {
        removed_F.row(i) = removedFacesList[i];
    }

    // Create new vertex list for MeshB_V excluding the vertices to remove
    MeshB_V.resize(MeshA_V.rows() - verticesToRemove.size(), MeshA_V.cols());
    int newIndex = 0;
    std::vector<int> oldToNewIndex(MeshA_V.rows(), -1);
    for (int v = 0; v < MeshA_V.rows(); ++v) {
        if (verticesToRemove.find(v) == verticesToRemove.end()) {
            MeshB_V.row(newIndex) = MeshA_V.row(v);
            oldToNewIndex[v] = newIndex++;
        }
    }

    // Create new face list for MeshB_F excluding faces that reference removed vertices
    std::vector<Eigen::Vector3i> newFaces;
    for (int f = 0; f < MeshA_F.rows(); ++f) {
        int v0 = MeshA_F(f, 0);
        int v1 = MeshA_F(f, 1);
        int v2 = MeshA_F(f, 2);

        if (verticesToRemove.find(v0) == verticesToRemove.end() &&
            verticesToRemove.find(v1) == verticesToRemove.end() &&
            verticesToRemove.find(v2) == verticesToRemove.end()) {
            newFaces.emplace_back(oldToNewIndex[v0], oldToNewIndex[v1], oldToNewIndex[v2]);
        }
    }
    MeshB_F = Eigen::MatrixXi(newFaces.size(), 3);
    for (size_t j = 0; j < newFaces.size(); ++j) {
        MeshB_F.row(j) = newFaces[j];
    }

    // Update removed_V with the coordinates of the removed vertices
    removed_V.resize(verticesToRemove.size(), 3);
    int index = 0;
    for (const auto& vertexIndex : verticesToRemove) {
        removed_V.row(index++) = MeshA_V.row(vertexIndex);
    }

    // Update border_V_new to include only the border vertices that were not removed
    std::vector<Eigen::RowVector3d> newBorderVertices;
    for (int i = 0; i < border_V.rows(); ++i) {
        bool isRemoved = false;
        for (int j = 0; j < removed_V.rows(); ++j) {
            if ((border_V.row(i) - removed_V.row(j)).norm() < tolerance) {
                isRemoved = true;
                break;
            }
        }
        if (!isRemoved) {
            newBorderVertices.push_back(border_V.row(i));
        }
    }
    border_V_new.resize(newBorderVertices.size(), 3);
    for (size_t j = 0; j < newBorderVertices.size(); ++j) {
        border_V_new.row(j) = newBorderVertices[j];
    }

    // Print output sizes
    std::cout << "Output MeshB_V size: " << MeshB_V.rows() << " vertices" << std::endl;
    std::cout << "Output MeshB_F size: " << MeshB_F.rows() << " faces" << std::endl;
    std::cout << "Output removed_V size: " << removed_V.rows() << " removed vertices" << std::endl;
    std::cout << "Output removed_F size: " << removed_F.rows() << " removed faces" << std::endl;
    std::cout << "Output border_V_new size: " << border_V_new.rows() << " border vertices (excluding removed)" << std::endl;
}


void showMeshAndPointCloud(const Eigen::MatrixXd& meshV, const Eigen::MatrixXi& meshF, const Eigen::MatrixXd& pointCloud) {
    // Initialize Polyscope
    polyscope::init();
    std::cout << "init Polyscope" << std::endl;

    // Register the mesh
    polyscope::registerSurfaceMesh("Mesh", meshV, meshF);

    // Register the point cloud
    polyscope::registerPointCloud("Point Cloud", pointCloud);

    // Show the Polyscope UI
    polyscope::show();
}

// Function to get the border vertices as an Eigen::MatrixXd
Eigen::MatrixXd getBorderVerticesMatrix(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    // Step 1: Compute the boundary loops using libigl's boundary_loop function
    std::vector<std::vector<int>> boundaryLoops;
    igl::boundary_loop(F, boundaryLoops);

    // Step 2: Collect all unique vertices that are part of any boundary loop
    std::unordered_set<int> uniqueBorderIndices;
    for (const auto& loop : boundaryLoops) {
        for (int vertex : loop) {
            uniqueBorderIndices.insert(vertex);
        }
    }

    // Step 3: Create the Eigen::MatrixXd for border vertices
    Eigen::MatrixXd borderVertices(uniqueBorderIndices.size(), V.cols());
    int row = 0;
    for (int index : uniqueBorderIndices) {
        borderVertices.row(row++) = V.row(index);
    }

    return borderVertices;
}


// Function to count the number of connected components in a mesh
int countConnectedComponents(const Eigen::MatrixXi& F) {
    // Step 1: Build adjacency list from the face matrix
    std::unordered_map<int, std::vector<int>> adjacencyList;
    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int v1 = F(i, j);
            int v2 = F(i, (j + 1) % 3);
            adjacencyList[v1].push_back(v2);
            adjacencyList[v2].push_back(v1);
        }
    }

    // Step 2: Perform DFS to find connected components
    std::unordered_map<int, bool> visited;
    for (const auto& [vertex, _] : adjacencyList) {
        visited[vertex] = false;
    }

    int numComponents = 0;
    for (const auto& [vertex, _] : adjacencyList) {
        if (!visited[vertex]) {
            // New component found, start DFS
            std::stack<int> stack;
            stack.push(vertex);
            visited[vertex] = true;

            while (!stack.empty()) {
                int v = stack.top();
                stack.pop();

                for (int neighbor : adjacencyList[v]) {
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        stack.push(neighbor);
                    }
                }
            }
            numComponents++;  // Finished traversing one component
        }
    }

    // Debug output
    std::cout << "Number of connected components: " << numComponents << std::endl;

    return numComponents;
}


void RemoveVerticesAndFaces(
    const Eigen::MatrixXd& Mesh1_V,   // Original vertices
    const Eigen::MatrixXi& Mesh1_F,   // Original faces
    const Eigen::MatrixXd& V3,        // Vertices to remove
    Eigen::MatrixXd& out_V,           // Output vertices after removal
    Eigen::MatrixXi& out_F,           // Output faces after removal
    Eigen::MatrixXd& removed_V,       // Removed vertices
    Eigen::MatrixXi& removed_F        // Removed faces
) {
    std::cout << "Input Vertices: " << Mesh1_V.rows() << std::endl;
    std::cout << "Input Faces: " << Mesh1_F.rows() << std::endl;
    std::cout << "Input Border Vertices: " << V3.rows() << std::endl;

    // Step 1: Identify indices of vertices to remove
    std::unordered_set<int> verticesToRemove;
    for (int i = 0; i < V3.rows(); ++i) {
        for (int j = 0; j < Mesh1_V.rows(); ++j) {
            if ((Mesh1_V.row(j) - V3.row(i)).norm() < 1e-6) {  // Threshold to match vertices
                verticesToRemove.insert(j);
                break;
            }
        }
    }

    // Step 2: Filter out faces that contain any of the vertices in `verticesToRemove`
    std::vector<Eigen::RowVector3i> newFaces;
    std::vector<Eigen::RowVector3i> removedFaces;
    for (int i = 0; i < Mesh1_F.rows(); ++i) {
        Eigen::RowVector3i face = Mesh1_F.row(i);
        if (verticesToRemove.count(face(0)) == 0 &&
            verticesToRemove.count(face(1)) == 0 &&
            verticesToRemove.count(face(2)) == 0) {
            newFaces.push_back(face);
        }
        else {
            removedFaces.push_back(face);
        }
    }

    // Step 3: Build new vertex list, reindex faces
    std::unordered_map<int, int> vertexMap; // Maps old index to new index
    std::vector<Eigen::RowVector3d> newVertices;
    std::vector<Eigen::RowVector3d> removedVertices;
    for (int i = 0; i < Mesh1_V.rows(); ++i) {
        if (verticesToRemove.count(i) == 0) {
            int newIndex = newVertices.size();
            newVertices.push_back(Mesh1_V.row(i));
            vertexMap[i] = newIndex;
        }
        else {
            removedVertices.push_back(Mesh1_V.row(i));
        }
    }

    // Update face indices to the new vertex indices
    for (auto& face : newFaces) {
        face(0) = vertexMap[face(0)];
        face(1) = vertexMap[face(1)];
        face(2) = vertexMap[face(2)];
    }

    // Step 4: Convert vectors to Eigen matrices for output
    out_V.resize(newVertices.size(), 3);
    for (int i = 0; i < newVertices.size(); ++i) {
        out_V.row(i) = newVertices[i];
    }

    out_F.resize(newFaces.size(), 3);
    for (int i = 0; i < newFaces.size(); ++i) {
        out_F.row(i) = newFaces[i];
    }

    removed_V.resize(removedVertices.size(), 3);
    for (int i = 0; i < removedVertices.size(); ++i) {
        removed_V.row(i) = removedVertices[i];
    }

    removed_F.resize(removedFaces.size(), 3);
    for (int i = 0; i < removedFaces.size(); ++i) {
        removed_F.row(i) = removedFaces[i];
    }

    std::cout << "Output Vertices: " << out_V.rows() << std::endl;
    std::cout << "Output Faces: " << out_F.rows() << std::endl;
    std::cout << "Removed Vertices: " << removed_V.rows() << std::endl;
    std::cout << "Removed Faces: " << removed_F.rows() << std::endl;
}














// Helper function to find the closest point on an edge (u, v) to a point p
Eigen::RowVector3d closestPointOnEdge(
    const Eigen::RowVector3d& u,
    const Eigen::RowVector3d& v,
    const Eigen::RowVector3d& p
) {
    Eigen::RowVector3d edge = v - u;
    double t = (p - u).dot(edge) / edge.squaredNorm();
    t = std::max(0.0, std::min(1.0, t));  // Clamp to [0, 1] to stay on the edge
    return u + t * edge;
}

// Helper function to extract a submesh given face indices
void extractSubmesh(
    const Eigen::MatrixXd& V,           // Original vertices
    const Eigen::MatrixXi& F,           // Original faces
    const std::vector<int>& faceIndices, // Indices of faces to include in submesh
    Eigen::MatrixXd& submesh_V,         // Output vertices of the submesh
    Eigen::MatrixXi& submesh_F          // Output faces of the submesh
) {
    std::unordered_map<int, int> vertexMap; // Maps original vertex indices to submesh indices
    std::vector<Eigen::RowVector3i> newFaces;
    std::vector<Eigen::RowVector3d> newVertices;

    // Iterate over the selected faces and build the new vertex and face lists
    for (int idx : faceIndices) {
        Eigen::RowVector3i face = F.row(idx);
        Eigen::RowVector3i newFace;

        // For each vertex in the face, map it to a new index in the submesh
        for (int j = 0; j < 3; ++j) {
            int originalIndex = face(j);
            if (vertexMap.count(originalIndex) == 0) {
                // If vertex not mapped yet, add it to the new vertex list
                int newIndex = newVertices.size();
                newVertices.push_back(V.row(originalIndex));
                vertexMap[originalIndex] = newIndex;
            }
            newFace(j) = vertexMap[originalIndex];
        }

        newFaces.push_back(newFace);
    }

    // Convert new vertices and faces to Eigen matrices
    submesh_V.resize(newVertices.size(), 3);
    for (int i = 0; i < newVertices.size(); ++i) {
        submesh_V.row(i) = newVertices[i];
    }

    submesh_F.resize(newFaces.size(), 3);
    for (int i = 0; i < newFaces.size(); ++i) {
        submesh_F.row(i) = newFaces[i];
    }

    std::cout << "Extracted submesh with " << submesh_V.rows() << " vertices and " << submesh_F.rows() << " faces.\n";
}

// Function to find connected components using DFS on the adjacency list
void findConnectedComponents(
    const Eigen::MatrixXi& F,
    std::vector<std::vector<int>>& components
) {
    // Step 1: Build the adjacency list
    std::unordered_map<int, std::vector<int>> adjacencyList;
    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int v1 = F(i, j);
            int v2 = F(i, (j + 1) % 3);
            adjacencyList[v1].push_back(v2);
            adjacencyList[v2].push_back(v1);
        }
    }

    // Step 2: Use DFS to find connected components
    std::unordered_map<int, bool> visited;
    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            visited[F(i, j)] = false;
        }
    }

    for (const auto& [vertex, _] : adjacencyList) {
        if (!visited[vertex]) {
            std::vector<int> component;
            std::stack<int> stack;
            stack.push(vertex);
            visited[vertex] = true;

            while (!stack.empty()) {
                int v = stack.top();
                stack.pop();
                component.push_back(v);

                for (int neighbor : adjacencyList[v]) {
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        stack.push(neighbor);
                    }
                }
            }

            components.push_back(component);
        }
    }

    std::cout << "Found " << components.size() << " connected components.\n";
}

//// Main function to split mesh by vertices on edges
//void splitMeshByVerticesOnEdges(
//    Eigen::MatrixXd& Mesh1_V,
//    Eigen::MatrixXi& Mesh1_F,
//    const Eigen::MatrixXd& V3,
//    Eigen::MatrixXd& MeshA_V,
//    Eigen::MatrixXi& MeshA_F,
//    Eigen::MatrixXd& MeshB_V,
//    Eigen::MatrixXi& MeshB_F
//) {
//    std::cout << "Starting mesh splitting...\n";
//    std::cout << "Original mesh has " << Mesh1_V.rows() << " vertices and " << Mesh1_F.rows() << " faces.\n";
//
//    // Step 1: Prepare containers for the new vertices and faces
//    Eigen::MatrixXd newVertices = Mesh1_V;
//    std::vector<Eigen::RowVector3i> newFaces;
//    std::set<int> splitVertices;
//
//    // Step 2: Iterate over each face in Mesh1_F to find edges that intersect with V3 vertices
//    for (int f = 0; f < Mesh1_F.rows(); ++f) {
//        Eigen::RowVector3i face = Mesh1_F.row(f);
//        for (int e = 0; e < 3; ++e) {
//            int v1 = face(e);
//            int v2 = face((e + 1) % 3);
//
//            Eigen::RowVector3d edgeStart = Mesh1_V.row(v1);
//            Eigen::RowVector3d edgeEnd = Mesh1_V.row(v2);
//
//            for (int v3Idx = 0; v3Idx < V3.rows(); ++v3Idx) {
//                Eigen::RowVector3d pointV3 = V3.row(v3Idx);
//                Eigen::RowVector3d closestPoint = closestPointOnEdge(edgeStart, edgeEnd, pointV3);
//                double distance = (closestPoint - pointV3).norm();
//
//                if (distance < 1e-6) {
//                    int newIndex = newVertices.rows();
//                    newVertices.conservativeResize(newVertices.rows() + 1, 3);
//                    newVertices.row(newIndex) = closestPoint;
//                    splitVertices.insert(newIndex);
//
//                    Eigen::RowVector3i newFace1, newFace2;
//                    newFace1 << v1, newIndex, face((e + 2) % 3);
//                    newFace2 << newIndex, v2, face((e + 2) % 3);
//                    newFaces.push_back(newFace1);
//                    newFaces.push_back(newFace2);
//
//                    std::cout << "Added split vertex at (" << closestPoint << ") between vertices " << v1 << " and " << v2 << ".\n";
//                    break;
//                }
//            }
//        }
//    }
//
//    // Step 3: Create new face matrix from updated face list
//    Mesh1_F.resize(newFaces.size(), 3);
//    for (int i = 0; i < newFaces.size(); ++i) {
//        Mesh1_F.row(i) = newFaces[i];
//    }
//    Mesh1_V = newVertices;
//    std::cout << "Updated mesh now has " << Mesh1_V.rows() << " vertices and " << Mesh1_F.rows() << " faces.\n";
//
//    // Step 4: Find connected components using DFS
//    std::vector<std::vector<int>> components;
//    findConnectedComponents(Mesh1_F, components);
//
//    // Step 5: Split into two submeshes based on component labels
//    if (components.size() >= 2) {
//        extractSubmesh(Mesh1_V, Mesh1_F, components[0], MeshA_V, MeshA_F);
//        extractSubmesh(Mesh1_V, Mesh1_F, components[1], MeshB_V, MeshB_F);
//        std::cout << "Mesh split into two components with " << MeshA_V.rows() << " and " << MeshB_V.rows() << " vertices.\n";
//    }
//    else {
//        std::cout << "Insufficient components found. Cannot split mesh into two parts.\n";
//    }
//}



// Function to split a mesh into two connected components
void splitMeshIn2(
    const Eigen::MatrixXd& V, // Input vertices
    const Eigen::MatrixXi& F, // Input faces
    Eigen::MatrixXd& V1,      // Output vertices for first component
    Eigen::MatrixXi& F1,      // Output faces for first component
    Eigen::MatrixXd& V2,      // Output vertices for second component
    Eigen::MatrixXi& F2       // Output faces for second component
) {
    // Step 1: Build adjacency list from faces
    std::vector<std::unordered_set<int>> adjacency(V.rows());
    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int v1 = F(i, j);
            int v2 = F(i, (j + 1) % 3);
            adjacency[v1].insert(v2);
            adjacency[v2].insert(v1);
        }
    }

    // Step 2: Initialize a visited array to mark which vertices have been processed
    std::vector<int> component(V.rows(), -1); // -1 indicates unvisited

    // Step 3: Helper function for BFS to collect a component
    auto bfs_component = [&](int start, int comp_id) -> std::vector<int> {
        std::vector<int> vertices;
        std::queue<int> queue;
        queue.push(start);
        component[start] = comp_id;

        while (!queue.empty()) {
            int v = queue.front();
            queue.pop();
            vertices.push_back(v);

            for (int neighbor : adjacency[v]) {
                if (component[neighbor] == -1) { // Not visited
                    component[neighbor] = comp_id;
                    queue.push(neighbor);
                }
            }
        }
        return vertices;
    };

    // Step 4: Find the two components
    std::vector<int> comp1, comp2;
    for (int i = 0; i < V.rows(); ++i) {
        if (component[i] == -1) {
            if (comp1.empty()) {
                comp1 = bfs_component(i, 0);
            }
            else {
                comp2 = bfs_component(i, 1);
                break; // We only need two components
            }
        }
    }

    // Step 5: Separate vertices and faces for each component
    std::unordered_map<int, int> map1, map2;
    std::vector<Eigen::RowVector3d> vertices1, vertices2;
    std::vector<Eigen::Vector3i> faces1, faces2;

    // Fill component 1 vertices and faces
    for (int v : comp1) {
        map1[v] = vertices1.size();
        vertices1.push_back(V.row(v));
    }
    for (int i = 0; i < F.rows(); ++i) {
        int v0 = F(i, 0), v1 = F(i, 1), v2 = F(i, 2);
        if (component[v0] == 0 && component[v1] == 0 && component[v2] == 0) {
            faces1.emplace_back(map1[v0], map1[v1], map1[v2]);
        }
    }

    // Fill component 2 vertices and faces
    for (int v : comp2) {
        map2[v] = vertices2.size();
        vertices2.push_back(V.row(v));
    }
    for (int i = 0; i < F.rows(); ++i) {
        int v0 = F(i, 0), v1 = F(i, 1), v2 = F(i, 2);
        if (component[v0] == 1 && component[v1] == 1 && component[v2] == 1) {
            faces2.emplace_back(map2[v0], map2[v1], map2[v2]);
        }
    }

    // Step 6: Convert to Eigen matrices
    V1.resize(vertices1.size(), 3);
    for (size_t i = 0; i < vertices1.size(); ++i) V1.row(i) = vertices1[i];

    F1.resize(faces1.size(), 3);
    for (size_t i = 0; i < faces1.size(); ++i) F1.row(i) = faces1[i];

    V2.resize(vertices2.size(), 3);
    for (size_t i = 0; i < vertices2.size(); ++i) V2.row(i) = vertices2[i];

    F2.resize(faces2.size(), 3);
    for (size_t i = 0; i < faces2.size(); ++i) F2.row(i) = faces2[i];

    std::cout << "Mesh successfully split into two components." << std::endl;
}


#include <iostream>
#include <set>
#include <utility>
#include <Eigen/Dense>

// Helper function to create a consistent representation of an edge
std::pair<int, int> createEdge(int v1, int v2) {
    return (v1 < v2) ? std::make_pair(v1, v2) : std::make_pair(v2, v1);
}

bool isVertexInRemoved(const Eigen::MatrixXd& removed_V, const Eigen::MatrixXd& mesh_V, int vertexIndex, double tolerance) {
    // Extract the vertex from mesh_V at the given vertexIndex
    Eigen::Vector3d vertex = mesh_V.row(vertexIndex);

    // Loop through all the vertices in removed_V to find a match
    for (int i = 0; i < removed_V.rows(); ++i) {
        // Extract the current vertex from removed_V
        Eigen::Vector3d removedVertex = removed_V.row(i);

        // Compare the difference for each coordinate, using the tolerance
        Eigen::Vector3d diff = vertex - removedVertex;

        // If all coordinates are within the tolerance, consider it a match
        if (std::abs(diff.x()) < tolerance &&
            std::abs(diff.y()) < tolerance &&
            std::abs(diff.z()) < tolerance) {
            return true;  // Vertex found in removed_V within tolerance
        }
    }

    // If no match was found, return false
    return false;
}


void FindMatchingEdges(
    const Eigen::MatrixXd& mesh_V,             // Vertex positions of the original mesh
    const Eigen::MatrixXi& mesh_F,             // Faces of the original mesh
    const Eigen::MatrixXd& removed_V,          // Vertex positions of the removed mesh
    const Eigen::MatrixXi& removed_F,          // Faces of the removed mesh
    double tolerance,                          // Tolerance for vertex matching
    Eigen::MatrixXd& mesh_V_restore,           // Vertex positions of the restored mesh
    Eigen::MatrixXi& mesh_F_restore            // Faces of the restored mesh
) {
    // Additional storage for added vertices and faces
    Eigen::MatrixXd added_V(0, 3);
    Eigen::MatrixXi added_F(0, 3);

    // Print input sizes
    std::cout << "Input Sizes:" << std::endl;
    std::cout << "mesh_V (original): " << mesh_V.rows() << " vertices" << std::endl;
    std::cout << "mesh_F (original): " << mesh_F.rows() << " faces" << std::endl;
    std::cout << "mesh_V_restore (initial): " << mesh_V_restore.rows() << " vertices" << std::endl;
    std::cout << "mesh_F_restore (initial): " << mesh_F_restore.rows() << " faces" << std::endl;
    std::cout << "removed_V: " << removed_V.rows() << " vertices" << std::endl;
    std::cout << "removed_F: " << removed_F.rows() << " faces" << std::endl;

    // Step 1: Build a set of edges from mesh_F
    std::set<std::pair<int, int>> meshEdges;
    for (int i = 0; i < mesh_F.rows(); ++i) {
        int v0 = mesh_F(i, 0);
        int v1 = mesh_F(i, 1);
        int v2 = mesh_F(i, 2);

        // Insert edges as ordered pairs
        meshEdges.insert(createEdge(v0, v1));
        meshEdges.insert(createEdge(v1, v2));
        meshEdges.insert(createEdge(v2, v0));
    }

    // Step 2: Find matching edges in removed_F and add missing vertices/faces
    for (int i = 0; i < removed_F.rows(); ++i) {
        int rv0 = removed_F(i, 0);
        int rv1 = removed_F(i, 1);
        int rv2 = removed_F(i, 2);

        bool foundMatch = false;

        auto findClosestVertex = [&](int vertexIndex) {
            for (int j = 0; j < mesh_V_restore.rows(); ++j) {


                Eigen::Vector3d diff = mesh_V.row(vertexIndex) - mesh_V_restore.row(j);


                // Compare individual coordinate differences instead of the Euclidean distance
                if (std::abs(diff.x()) < tolerance &&
                    std::abs(diff.y()) < tolerance &&
                    std::abs(diff.z()) < tolerance) {
                    std::cout << "compared and accepted : " << mesh_V.row(vertexIndex) << " to " << mesh_V_restore.row(j)  << 
                        "at j: " << j << std::endl;
                    return j;
                }
            }
            return -1;
        };

        // Check and add faces for each edge
        auto processEdge = [&](int idx0, int idx1, int idx2, int rv2) {
            if (idx0 != -1 && idx1 != -1) {
                if (idx2 == -1) {
                    // Check if the vertex is in removed_V before adding
                    if (isVertexInRemoved(removed_V, mesh_V, rv2, tolerance) || true) {
                        // Add missing vertex
                        mesh_V_restore.conservativeResize(mesh_V_restore.rows() + 1, Eigen::NoChange);
                        mesh_V_restore.row(mesh_V_restore.rows() - 1) = mesh_V.row(rv2);
                        idx2 = mesh_V_restore.rows() - 1;
                        added_V.conservativeResize(added_V.rows() + 1, Eigen::NoChange);
                        added_V.row(added_V.rows() - 1) = mesh_V.row(rv2);
                        std::cout << "Added vertex " << idx2 << " to mesh_V_restore: " << mesh_V_restore.row(idx2).transpose() << std::endl;
                    }
                    else {
                        // Vertex not found in removed_V, skip adding face
                        std::cout << "Skipped adding vertex " << rv2 << " and corresponding face as it is not in removed_V." << std::endl;
                        return;
                    }
                }

                // Add new face
                mesh_F_restore.conservativeResize(mesh_F_restore.rows() + 1, Eigen::NoChange);
                mesh_F_restore.row(mesh_F_restore.rows() - 1) << idx0, idx1, idx2;
                added_F.conservativeResize(added_F.rows() + 1, Eigen::NoChange);
                added_F.row(added_F.rows() - 1) << idx0, idx1, idx2;
                std::cout << "Added face [" << idx0 << ", " << idx1 << ", " << idx2 << "] to mesh_F_restore." << std::endl;
            }
        };

        // Process each edge
        if (meshEdges.find(createEdge(rv0, rv1)) != meshEdges.end()) {
            std::cout << "1||| rv0, rv1, rv2 " << rv0 << ", " << rv1 << ", " << rv2 << std::endl;
            int idx0 = findClosestVertex(rv0);
            int idx1 = findClosestVertex(rv1);
            int idx2 = findClosestVertex(rv2);
            processEdge(idx0, idx1, idx2, rv2);
            foundMatch = true;
        }

        if (meshEdges.find(createEdge(rv1, rv2)) != meshEdges.end()) {
            std::cout << "2 |||| rv0, rv1, rv2 " << rv0 << ", " << rv1 << ", " << rv2 << std::endl;

            int idx1 = findClosestVertex(rv1);
            int idx2 = findClosestVertex(rv2);
            int idx0 = findClosestVertex(rv0);
            processEdge(idx1, idx2, idx0, rv0);
            foundMatch = true;
        }

        if (meshEdges.find(createEdge(rv2, rv0)) != meshEdges.end()) {
            std::cout << "3 ||| rv0, rv1, rv2 " << rv0 << ", " << rv1 << ", " << rv2 << std::endl;

            int idx2 = findClosestVertex(rv2);
            int idx0 = findClosestVertex(rv0);
            int idx1 = findClosestVertex(rv1);
            processEdge(idx2, idx0, idx1, rv1);
            foundMatch = true;
        }

        if (!foundMatch) {
            std::cout << "No matching edge found for removed face " << i << std::endl;
        }
    }

    // Print final results
    std::cout << "OUTPUT Sizes:" << std::endl;
    std::cout << "mesh_V_restore (final): " << mesh_V_restore.rows() << " vertices" << std::endl;
    std::cout << "mesh_F_restore (final): " << mesh_F_restore.rows() << " faces" << std::endl;
    std::cout << "Vertices added during restoration:\n" << added_V << std::endl;
    std::cout << "Vertices added during restoration (length):" << added_V.rows() << std::endl;
    std::cout << "Faces added during restoration: (length)" << added_F.rows() << std::endl;
    std::cout << "Faces added during restoration:\n" << added_F << std::endl;
    std::cout << "Restoration complete." << std::endl;
}




