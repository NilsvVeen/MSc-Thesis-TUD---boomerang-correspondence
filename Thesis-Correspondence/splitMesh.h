#ifndef SPLIT_MESH_H
#define SPLIT_MESH_H

#include <Eigen/Dense>
#include <set>
#include <vector>
#include <unordered_map>
#include <stack>

// Function to find the closest point on an edge to a given point
Eigen::RowVector3d closestPointOnEdge(
    const Eigen::RowVector3d& u,
    const Eigen::RowVector3d& v,
    const Eigen::RowVector3d& p
);

// Function to extract a submesh given face indices
void extractSubmesh(
    const Eigen::MatrixXd& V,           // Original vertices
    const Eigen::MatrixXi& F,           // Original faces
    const std::vector<int>& faceIndices, // Indices of faces to include in submesh
    Eigen::MatrixXd& submesh_V,         // Output vertices of the submesh
    Eigen::MatrixXi& submesh_F          // Output faces of the submesh
);

// Function to find connected components using DFS
void findConnectedComponents(
    const Eigen::MatrixXi& F,
    std::vector<std::vector<int>>& components
);

// Main function to split mesh by vertices on edges
void splitMeshByVerticesOnEdges(
    Eigen::MatrixXd& Mesh1_V,
    Eigen::MatrixXi& Mesh1_F,
    const Eigen::MatrixXd& V3,
    Eigen::MatrixXd& MeshA_V,
    Eigen::MatrixXi& MeshA_F,
    Eigen::MatrixXd& MeshB_V,
    Eigen::MatrixXi& MeshB_F
);

void RemoveVerticesAndFaces(
    const Eigen::MatrixXd& Mesh1_V,   // Original vertices
    const Eigen::MatrixXi& Mesh1_F,   // Original faces
    const Eigen::MatrixXd& V3,        // Vertices to remove
    Eigen::MatrixXd& out_V,           // Output vertices after removal
    Eigen::MatrixXi& out_F            // Output faces after removal
);


int countConnectedComponents(const Eigen::MatrixXi& F);
Eigen::MatrixXd getBorderVerticesMatrix(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

void showMeshAndPointCloud(const Eigen::MatrixXd& meshV, const Eigen::MatrixXi& meshF, const Eigen::MatrixXd& pointCloud);

void removeVerticesWithTwoFacesAndBorderEdges(
    const Eigen::MatrixXd& MeshA_V,
    const Eigen::MatrixXi& MeshA_F,
    Eigen::MatrixXd& border_V, // Change border_V to Eigen::MatrixXd
    Eigen::MatrixXd& MeshB_V,
    Eigen::MatrixXi& MeshB_F,
    Eigen::MatrixXd& removed_V // Changed to Eigen::MatrixXd
);


#endif // SPLIT_MESH_H
