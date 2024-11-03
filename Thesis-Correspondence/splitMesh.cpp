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

// Function to remove vertices with specific conditions
void removeVerticesWithTwoFacesAndBorderEdges(
    const Eigen::MatrixXd& MeshA_V,
    const Eigen::MatrixXi& MeshA_F,
    Eigen::MatrixXd& border_V, // Matrix of border vertices
    Eigen::MatrixXd& MeshB_V,
    Eigen::MatrixXi& MeshB_F,
    Eigen::MatrixXd& removed_V // Matrix of removed vertices
) {

    // Print input sizes
    std::cout << "Input MeshA_V size: " << MeshA_V.rows() << " vertices" << std::endl;
    std::cout << "Input MeshA_F size: " << MeshA_F.rows() << " faces" << std::endl;
    std::cout << "Input border_V size: " << border_V.rows() << " border vertices" << std::endl;


    // Create a map to count the number of faces each vertex belongs to
    std::unordered_map<int, int> vertexFaceCount;

    // Count the faces for each vertex in MeshA_F
    for (int f = 0; f < MeshA_F.rows(); ++f) {
        for (int j = 0; j < 3; ++j) {
            vertexFaceCount[MeshA_F(f, j)]++;
        }
    }

    // Use a set to keep track of vertices to remove
    std::unordered_set<int> verticesToRemove;
    // Use a set to keep track of printed border vertices
    std::unordered_set<int> printedVertices;

    // Loop over border vertices and check their face count
    for (int i = 0; i < border_V.rows(); ++i) {
        int borderVertexIndex = static_cast<int>(border_V(i, 0)); // Assuming the first column has the vertex indices

        // Check if the border vertex is part of exactly 2 faces
        if (vertexFaceCount[borderVertexIndex] == 2) {
            // Check if this vertex has already been printed
            if (printedVertices.find(borderVertexIndex) == printedVertices.end()) {
                // Retrieve the vertex coordinates from MeshA_V
                Eigen::RowVector3d vertexCoordinates = MeshA_V.row(borderVertexIndex);

                // Print the border vertex index and its coordinates
                std::cout << "Border vertex " << borderVertexIndex
                    << " (MeshA_V index) at coordinates "
                    << vertexCoordinates.transpose()
                    << " is part of 2 faces." << std::endl;

                printedVertices.insert(borderVertexIndex); // Add to printed set
                verticesToRemove.insert(borderVertexIndex); // Mark this vertex for removal
            }
        }
    }

    // Create new vertex list for MeshB_V excluding the vertices to remove
    std::vector<Eigen::RowVector3d> newVertices;
    for (int v = 0; v < MeshA_V.rows(); ++v) {
        if (verticesToRemove.find(v) == verticesToRemove.end()) {
            newVertices.push_back(MeshA_V.row(v));
        }
    }
    MeshB_V = Eigen::MatrixXd(newVertices.size(), 3);
    for (size_t j = 0; j < newVertices.size(); ++j) {
        MeshB_V.row(j) = newVertices[j];
    }

    // Create new face list for MeshB_F excluding faces that reference removed vertices
    std::vector<Eigen::Vector3i> newFaces;
    for (int f = 0; f < MeshA_F.rows(); ++f) {
        bool faceValid = true;
        for (int j = 0; j < 3; ++j) {
            if (verticesToRemove.find(MeshA_F(f, j)) != verticesToRemove.end()) {
                faceValid = false;
                break;
            }
        }
        if (faceValid) {
            newFaces.push_back(MeshA_F.row(f));
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
        removed_V.row(index) = MeshA_V.row(vertexIndex);
        index++;
    }

    // Update border_V to exclude the removed vertices
    std::vector<Eigen::RowVector3d> newBorderVertices;
    for (int i = 0; i < border_V.rows(); ++i) {
        int borderVertexIndex = static_cast<int>(border_V(i, 0));
        if (verticesToRemove.find(borderVertexIndex) == verticesToRemove.end()) {
            newBorderVertices.push_back(border_V.row(i));
        }
    }
    border_V.resize(newBorderVertices.size(), 3);
    for (size_t j = 0; j < newBorderVertices.size(); ++j) {
        border_V.row(j) = newBorderVertices[j];
    }

    // Print output sizes
    std::cout << "Output MeshB_V size: " << MeshB_V.rows() << " vertices" << std::endl;
    std::cout << "Output MeshB_F size: " << MeshB_F.rows() << " faces" << std::endl;
    std::cout << "Output removed_V size: " << removed_V.rows() << " removed vertices" << std::endl;
    std::cout << "Output border_V size: " << border_V.rows() << " border vertices" << std::endl;
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
    Eigen::MatrixXi& out_F            // Output faces after removal
) {

    std::cout << "input Vertices: " << Mesh1_V.rows() << std::endl;
    std::cout << "input Faces: " << Mesh1_F.rows() << std::endl;
    std::cout << "input border vertices: " << V3.rows() << std::endl;

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
    for (int i = 0; i < Mesh1_F.rows(); ++i) {
        Eigen::RowVector3i face = Mesh1_F.row(i);
        if (verticesToRemove.count(face(0)) == 0 &&
            verticesToRemove.count(face(1)) == 0 &&
            verticesToRemove.count(face(2)) == 0) {
            newFaces.push_back(face);
        }
    }

    // Step 3: Build new vertex list, reindex faces
    std::unordered_map<int, int> vertexMap; // Maps old index to new index
    std::vector<Eigen::RowVector3d> newVertices;
    for (int i = 0; i < Mesh1_V.rows(); ++i) {
        if (verticesToRemove.count(i) == 0) {
            int newIndex = newVertices.size();
            newVertices.push_back(Mesh1_V.row(i));
            vertexMap[i] = newIndex;
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


    std::cout << "output Vertices: " << out_V.rows() << std::endl;
    std::cout << "output Faces: " << out_F.rows() << std::endl;
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

// Main function to split mesh by vertices on edges
void splitMeshByVerticesOnEdges(
    Eigen::MatrixXd& Mesh1_V,
    Eigen::MatrixXi& Mesh1_F,
    const Eigen::MatrixXd& V3,
    Eigen::MatrixXd& MeshA_V,
    Eigen::MatrixXi& MeshA_F,
    Eigen::MatrixXd& MeshB_V,
    Eigen::MatrixXi& MeshB_F
) {
    std::cout << "Starting mesh splitting...\n";
    std::cout << "Original mesh has " << Mesh1_V.rows() << " vertices and " << Mesh1_F.rows() << " faces.\n";

    // Step 1: Prepare containers for the new vertices and faces
    Eigen::MatrixXd newVertices = Mesh1_V;
    std::vector<Eigen::RowVector3i> newFaces;
    std::set<int> splitVertices;

    // Step 2: Iterate over each face in Mesh1_F to find edges that intersect with V3 vertices
    for (int f = 0; f < Mesh1_F.rows(); ++f) {
        Eigen::RowVector3i face = Mesh1_F.row(f);
        for (int e = 0; e < 3; ++e) {
            int v1 = face(e);
            int v2 = face((e + 1) % 3);

            Eigen::RowVector3d edgeStart = Mesh1_V.row(v1);
            Eigen::RowVector3d edgeEnd = Mesh1_V.row(v2);

            for (int v3Idx = 0; v3Idx < V3.rows(); ++v3Idx) {
                Eigen::RowVector3d pointV3 = V3.row(v3Idx);
                Eigen::RowVector3d closestPoint = closestPointOnEdge(edgeStart, edgeEnd, pointV3);
                double distance = (closestPoint - pointV3).norm();

                if (distance < 1e-6) {
                    int newIndex = newVertices.rows();
                    newVertices.conservativeResize(newVertices.rows() + 1, 3);
                    newVertices.row(newIndex) = closestPoint;
                    splitVertices.insert(newIndex);

                    Eigen::RowVector3i newFace1, newFace2;
                    newFace1 << v1, newIndex, face((e + 2) % 3);
                    newFace2 << newIndex, v2, face((e + 2) % 3);
                    newFaces.push_back(newFace1);
                    newFaces.push_back(newFace2);

                    std::cout << "Added split vertex at (" << closestPoint << ") between vertices " << v1 << " and " << v2 << ".\n";
                    break;
                }
            }
        }
    }

    // Step 3: Create new face matrix from updated face list
    Mesh1_F.resize(newFaces.size(), 3);
    for (int i = 0; i < newFaces.size(); ++i) {
        Mesh1_F.row(i) = newFaces[i];
    }
    Mesh1_V = newVertices;
    std::cout << "Updated mesh now has " << Mesh1_V.rows() << " vertices and " << Mesh1_F.rows() << " faces.\n";

    // Step 4: Find connected components using DFS
    std::vector<std::vector<int>> components;
    findConnectedComponents(Mesh1_F, components);

    // Step 5: Split into two submeshes based on component labels
    if (components.size() >= 2) {
        extractSubmesh(Mesh1_V, Mesh1_F, components[0], MeshA_V, MeshA_F);
        extractSubmesh(Mesh1_V, Mesh1_F, components[1], MeshB_V, MeshB_F);
        std::cout << "Mesh split into two components with " << MeshA_V.rows() << " and " << MeshB_V.rows() << " vertices.\n";
    }
    else {
        std::cout << "Insufficient components found. Cannot split mesh into two parts.\n";
    }
}
