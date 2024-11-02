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
    Eigen::MatrixXd& border_V, // Change border_V to Eigen::MatrixXd
    Eigen::MatrixXd& MeshB_V,
    Eigen::MatrixXi& MeshB_F,
    Eigen::MatrixXd& removed_V // Changed to Eigen::MatrixXd
) {

    std::cout << "border vertices old: " << border_V.rows() << std::endl;


    // Step 1: Create an adjacency list to keep track of vertex connectivity
    std::vector<std::vector<int>> vertexFaces(MeshA_V.rows());
    for (int f = 0; f < MeshA_F.rows(); ++f) {
        for (int j = 0; j < 3; ++j) {
            vertexFaces[MeshA_F(f, j)].push_back(f);
        }
    }

    // Convert border_V to a set for fast lookups
    std::unordered_set<int> borderSet;
    for (int i = 0; i < border_V.rows(); ++i) {
        borderSet.insert(static_cast<int>(border_V(i, 0))); // Assuming border_V contains vertex indices in the first column
    }

    // Step 2: Find vertices to remove
    std::unordered_set<int> toRemove;
    for (int v = 0; v < MeshA_V.rows(); ++v) {
        if (vertexFaces[v].size() == 2) {
            int face1 = vertexFaces[v][0];
            int face2 = vertexFaces[v][1];

            // Get the other vertices of the two faces
            std::unordered_set<int> neighboringVertices;
            for (int j = 0; j < 3; ++j) {
                if (MeshA_F(face1, j) != v) neighboringVertices.insert(MeshA_F(face1, j));
                if (MeshA_F(face2, j) != v) neighboringVertices.insert(MeshA_F(face2, j));
            }

            // Check if both neighboring vertices are in the border set
            if (neighboringVertices.size() == 2) {
                bool bothInBorder = true;
                for (int nv : neighboringVertices) {
                    if (borderSet.find(nv) == borderSet.end()) {
                        bothInBorder = false;
                        break;
                    }
                }
                if (bothInBorder) {
                    toRemove.insert(v);
                }
            }
        }
    }

    // Step 3: Construct new mesh and list of removed vertices
    std::unordered_map<int, int> oldToNewIndex; // Map old vertex indices to new indices
    std::vector<Eigen::RowVector3d> removedVertices; // To collect removed vertices
    for (int v = 0; v < MeshA_V.rows(); ++v) {
        if (toRemove.find(v) == toRemove.end()) { // Not marked for removal
            int newIndex = MeshB_V.rows();
            MeshB_V.conservativeResize(newIndex + 1, 3);
            MeshB_V.row(newIndex) = MeshA_V.row(v);
            oldToNewIndex[v] = newIndex;
        }
        else {
            // Collect the removed vertex
            removedVertices.push_back(MeshA_V.row(v));
        }
    }

    // Convert to Eigen matrix for removed vertices
    removed_V.resize(removedVertices.size(), 3);
    for (size_t i = 0; i < removedVertices.size(); ++i) {
        removed_V.row(i) = removedVertices[i];
    }

    // Step 4: Create new faces for MeshB
    std::vector<Eigen::RowVector3i> newFaces;
    for (int f = 0; f < MeshA_F.rows(); ++f) {
        if (toRemove.find(MeshA_F(f, 0)) == toRemove.end() &&
            toRemove.find(MeshA_F(f, 1)) == toRemove.end() &&
            toRemove.find(MeshA_F(f, 2)) == toRemove.end()) {
            Eigen::RowVector3i newFace;
            for (int j = 0; j < 3; ++j) {
                newFace(j) = oldToNewIndex[MeshA_F(f, j)];
            }
            newFaces.push_back(newFace);
        }
    }

    // Convert to Eigen matrix for faces
    MeshB_F.resize(newFaces.size(), 3);
    for (size_t i = 0; i < newFaces.size(); ++i) {
        MeshB_F.row(i) = newFaces[i];
    }

    // Step 5: Remove the vertices from border_V
    std::vector<int> newBorderVertices;
    for (int i = 0; i < border_V.rows(); ++i) {
        if (toRemove.find(static_cast<int>(border_V(i, 0))) == toRemove.end()) {
            newBorderVertices.push_back(static_cast<int>(border_V(i, 0)));
        }
    }

    // Resize border_V to contain only the new border vertices
    border_V.resize(newBorderVertices.size(), 1);
    for (size_t i = 0; i < newBorderVertices.size(); ++i) {
        border_V(i, 0) = newBorderVertices[i];
    }

    std::cout << "removed border vertices: " << removed_V.rows() << std::endl;
    std::cout << "border vertices new: " << border_V.rows() << std::endl;
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
