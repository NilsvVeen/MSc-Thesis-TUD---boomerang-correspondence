#include <igl/readOFF.h>
#include <igl/lscm.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <iostream>
#include "stl_utils.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <regex>
#include <Eigen/Dense>
#include <algorithm>
#include <igl/harmonic.h>
#include <igl/edges.h>
#include <igl/flip_avoiding_line_search.h>
#include <igl/arap.h>

// Function to read and concatenate point clouds in numeric order from files matching a pattern
Eigen::MatrixXd readAndConcatenatePointClouds(const std::string& folderPath, const std::string& filePatternStr) {
    std::vector<std::pair<int, std::string>> orderedFiles;
    std::regex filePattern(filePatternStr);

    // Collect and parse filenames
    for (const auto& entry : std::filesystem::directory_iterator(folderPath)) {
        if (entry.is_regular_file()) {
            std::smatch match;
            std::string filename = entry.path().filename().string();
            if (std::regex_match(filename, match, filePattern)) {
                int fileIndex = std::stoi(match[1].str()); // Extract numeric index
                orderedFiles.emplace_back(fileIndex, entry.path().string());
                std::cout << "Found matching file: " << filename << " with index " << fileIndex << std::endl;
            }
        }
    }

    // Sort files by extracted numeric index
    std::sort(orderedFiles.begin(), orderedFiles.end());
    std::cout << "Files sorted in order of index." << std::endl;

    // Concatenate vertices from each file
    std::vector<Eigen::MatrixXd> matrices;
    int totalRows = 0;

    for (const auto& [index, filepath] : orderedFiles) {
        std::cout << "Reading vertices from file: " << filepath << std::endl;
        Eigen::MatrixXd V = readVerticesFromPLY(filepath);
        if (V.size() > 0) {
            totalRows += V.rows();
            matrices.push_back(V);
            std::cout << "File " << filepath << " has " << V.rows() << " vertices." << std::endl;
        }
        else {
            std::cerr << "Warning: No vertices read from file " << filepath << std::endl;
        }
    }

    // Combine all matrices into one
    Eigen::MatrixXd V3(totalRows, 3);
    int currentRow = 0;
    for (const auto& mat : matrices) {
        V3.block(currentRow, 0, mat.rows(), mat.cols()) = mat;
        std::cout << "Copying " << mat.rows() << " rows into V3 at position " << currentRow << std::endl;
        currentRow += mat.rows();
    }

    std::cout << "Total rows in concatenated matrix V3: " << V3.rows() << std::endl;
    return V3;
}

void CalculateGenus(Eigen::MatrixXd V, Eigen::MatrixXi F) {
    // Count vertices and faces
    int num_vertices = V.rows();
    int num_faces = F.rows();

    // Compute edges
    Eigen::MatrixXi E; // Edges
    igl::edges(F, E);
    int num_edges = E.rows();

    // Calculate Euler characteristic
    int euler_characteristic = num_vertices - num_edges + num_faces;

    // Calculate genus
    int genus = 1 - (euler_characteristic / 2);

    // Output results
    std::cout << "Vertices: " << num_vertices << std::endl;
    std::cout << "Edges: " << num_edges << std::endl;
    std::cout << "Faces: " << num_faces << std::endl;
    std::cout << "Euler Characteristic: " << euler_characteristic << std::endl;
    std::cout << "Genus: " << genus << std::endl;
}


//// Function to compute the boundary smoothness energy with stronger smoothing
//Eigen::VectorXd computeBoundaryEnergy(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXi& B, const Eigen::MatrixXd& UV) {
//    Eigen::VectorXd energy(B.size());
//
//    // Compute the smoothness energy based on neighboring boundary UVs
//    for (int i = 0; i < B.size(); ++i) {
//        int vIdx = B[i];
//        Eigen::Vector2d uv = UV.row(vIdx);  // This is a 1x2 row vector
//
//        // Neighbors: Wrap around using modulus to avoid index out of bounds
//        Eigen::Vector2d neighborUV = UV.row(B[(i + 1) % B.size()]).transpose();  // Next neighbor
//        Eigen::Vector2d prevUV = UV.row(B[(i - 1 + B.size()) % B.size()]).transpose();  // Previous neighbor
//
//        // Compute the smoothness energy as the difference from the previous and next boundary UVs
//        Eigen::Vector2d deltaNext = uv - neighborUV;
//        Eigen::Vector2d deltaPrev = uv - prevUV;
//
//        // Use the combined difference as an energy measure
//        energy(i) = (deltaNext.norm() + deltaPrev.norm()) / 2.0;  // Average energy based on both neighbors
//    }
//
//    return energy;
//}
//
//// Modified LSCM with stronger smoothness energy on the boundary
//void lscmWithSmoothBoundary(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXi& B, const Eigen::MatrixXd& BC, Eigen::MatrixXd& UV) {
//    // Perform LSCM first
//    igl::lscm(V, F, B, BC, UV);
//
//    // After LSCM, add smoothness energy adjustment
//    Eigen::VectorXd boundaryEnergy = computeBoundaryEnergy(V, F, B, UV);
//
//    // Stronger smoothing iteration to reduce large peaks
//    int numIterations = 200;  // Number of iterations to smooth the boundary
//    double smoothingFactor = 1;  // Strength of smoothing
//
//    for (int iteration = 0; iteration < numIterations; ++iteration) {
//        for (int i = 0; i < B.size(); ++i) {
//            int vIdx = B[i];
//
//            // Neighbors: Wrap around using modulus to avoid index out of bounds
//            Eigen::Vector2d neighborUV = UV.row(B[(i + 1) % B.size()]).transpose();  // Next neighbor
//            Eigen::Vector2d prevUV = UV.row(B[(i - 1 + B.size()) % B.size()]).transpose();  // Previous neighbor
//
//            // Laplacian smoothing: Take the average of the neighbors
//            Eigen::Vector2d smoothedUV = (prevUV + neighborUV) / 2.0;
//
//            // Apply smoothing with a weighted factor
//            UV.row(vIdx) = (1 - smoothingFactor) * UV.row(vIdx) + smoothingFactor * smoothedUV.transpose();
//        }
//    }
//
//    // Optionally, re-run LSCM after adjustment to further refine the result
//    igl::lscm(V, F, B, BC, UV);
//}



bool paramsurface5(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& UV, const Eigen::MatrixXd& boundary_vertices, bool boundaryEnabled, const Eigen::MatrixXd& boundary_vertices_other)
{
    // Print size of V and F
    std::cout << "Input V size: " << V.rows() << " x " << V.cols() << std::endl;
    std::cout << "Input F size: " << F.rows() << " x " << F.cols() << std::endl;

    //// Print the first 5 rows of V (if there are at least 5)
    //std::cout << "First 5 rows of V:" << std::endl;
    //for (int i = 0; i < std::min(5, static_cast<int>(V.rows())); ++i) {
    //    std::cout << "V(" << i << ") = [";
    //    for (int j = 0; j < V.cols(); ++j) {
    //        std::cout << V(i, j);
    //        if (j < V.cols() - 1) std::cout << ", ";
    //    }
    //    std::cout << "]" << std::endl;
    //}

    //// Print the first 5 rows of F (if there are at least 5)
    //std::cout << "First 5 rows of F:" << std::endl;
    //for (int i = 0; i < std::min(5, static_cast<int>(F.rows())); ++i) {
    //    std::cout << "F(" << i << ") = [";
    //    for (int j = 0; j < F.cols(); ++j) {
    //        std::cout << F(i, j);
    //        if (j < F.cols() - 1) std::cout << ", ";
    //    }
    //    std::cout << "]" << std::endl;
    //}

    UV.resize(V.rows(), 2);  // UV coordinates need two dimensions (u, v)

    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < F.cols(); ++j) {
            if (F(i, j) < 0 || F(i, j) >= V.rows()) {
                std::cout << "F(" << i << ", " << j << ") = " << F(i, j) << std::endl;
                std::cerr << "Error: Invalid face index at F(" << i << ", " << j << ")." << std::endl;
                return false;
            }
        }
    }

    std::cout << "vertices input for boudnary " << boundary_vertices.rows() << std::endl;


    // Initialize Polyscope
    polyscope::init();

    // Register the 3D mesh
    polyscope::registerSurfaceMesh("3D Mesh before", V, F);

    polyscope::show();

    CalculateGenus(V, F);

    // Perform LSCM parametrization
    std::cout << "Performing LSCM parametrization..." << std::endl;

    Eigen::MatrixXd boundaryUV(boundary_vertices.rows(), 2);

    if (boundaryEnabled) {
        // Initialize boundary constraints
        Eigen::MatrixXd BC(boundary_vertices.rows(), 2);
        Eigen::VectorXi B(boundary_vertices.rows());

        for (int i = 0; i < boundary_vertices.rows(); ++i) {
            double min_distance = std::numeric_limits<double>::max();
            int closest_vertex = -1;

            // Find the closest vertex in V to the current boundary vertex
            for (int j = 0; j < V.rows(); ++j) {
                double distance = (V.row(j) - boundary_vertices.row(i)).squaredNorm();
                if (distance < min_distance) {
                    min_distance = distance;
                    closest_vertex = j;
                }
            }

            if (closest_vertex == -1) {
                std::cerr << "Error: No closest vertex found for boundary vertex at index " << i
                    << " (" << boundary_vertices.row(i) << ")" << std::endl;
                return false; // Early exit if no closest vertex was found
            }

            B(i) = closest_vertex;  // Assign the closest vertex index

            // Set UV coordinates in BC
            //bool CircleUV = false;  // Toggle between different UV generation methods
            //if (CircleUV) {
            //    double angle = 2.0 * 3.14159265358979323846  * i / boundary_vertices.rows();
            //    BC(i, 0) = std::cos(angle);  // u-coordinate
            //    BC(i, 1) = std::sin(angle);  // v-coordinate
            //}

            bool ProjectionUV = true;  // If ProjectionUV is enabled
            if (ProjectionUV) {
                BC(i, 0) = boundary_vertices_other(i, 0);  // u-coordinate
                BC(i, 1) = boundary_vertices_other(i, 1);  // v-coordinate
                boundaryUV(i, 0) = BC(i, 0);
                boundaryUV(i, 1) = BC(i, 1);
            }
        }

        // Debug: Visualize the boundary vertices (B) and their UV constraints (BC)
        for (int i = 0; i < B.rows(); ++i) {
            std::cout << "Boundary Vertex Index: " << B(i)
                << ", UV: (" << BC(i, 0) << ", " << BC(i, 1) << ")" << std::endl;
        }


        // Perform LSCM with additional parameters
        igl::lscm(V, F, B, BC, UV);
        //igl::harmonic(V, F, B, BC, 2, UV); // The last parameter is the harmonic order
        //lscmWithSmoothBoundary(V, F, B, BC, UV);




    }
    else {
        igl::lscm(V, F, UV);
    }







    CalculateGenus(V, F);


    // Print size of UV matrix
    std::cout << "UV matrix size after LSCM: " << UV.rows() << " x " << UV.cols() << std::endl;

    // Check if the UV matrix is valid (it shouldn't be empty)
    if (UV.size() == 0)
    {
        std::cerr << "Failed to compute UV parametrization!" << std::endl;
        return false;
    }

    std::cout << "UV parametrization computed successfully." << std::endl;

    // Add a small constant z value to the UV coordinates to create a 3D version
    Eigen::MatrixXd UV3D(UV.rows(), 3);
    double z_value = 0.1; // Adjust this value for how far above the UV plane you want it
    for (int i = 0; i < UV.rows(); ++i) {
        UV3D(i, 0) = UV(i, 0) ; // u coordinate
        UV3D(i, 1) = UV(i, 1); // v coordinate
        //UV3D(i, 0) = UV(i, 0) * 100; // u coordinate
        //UV3D(i, 1) = UV(i, 1) * 100; // v coordinate
        UV3D(i, 2) = z_value;  // constant z value
    }

    writeVerticesToPLY("aaaaaaaaaaaaaaaaaaaaaaaaaaaaa", UV3D);


    // Initialize Polyscope
    polyscope::init();

    // Register the 3D mesh
    auto* mesh3D = polyscope::registerSurfaceMesh("3D Mesh", V, F);

    // Register the UV-mapped 3D mesh (with z-coordinate)
    auto* meshUV = polyscope::registerSurfaceMesh("UV Mapped Mesh", UV3D, F);

    // Add vertex parameterization to the 3D mesh
    mesh3D->addVertexParameterizationQuantity("LSCM parameterization", UV);

    // Set some options for the 3D mesh visualization
    mesh3D->setEnabled(true);  // Start with the 3D mesh enabled
    mesh3D->setSmoothShade(true);
    mesh3D->setEdgeWidth(1.0); // Optional: control edge display for 3D mesh

    // Set options for the UV mesh visualization
    meshUV->setEnabled(false);  // UV mesh is disabled initially
    meshUV->setEdgeWidth(1.0);  // Adjust edge width for the UV mesh visualization

    // Add a callback to toggle between the 3D and UV mesh views
    polyscope::state::userCallback = [&]()
    {
        if (ImGui::Button("Show 3D Mesh"))
        {
            mesh3D->setEnabled(true);
            meshUV->setEnabled(false);
        }
        if (ImGui::Button("Show UV Mapped Mesh"))
        {
            mesh3D->setEnabled(false);
            meshUV->setEnabled(true);
        }
    };

    auto* meshUVaaa = polyscope::registerPointCloud2D("boundary", boundaryUV);

    CalculateGenus(V, F);













    // Launch Polyscope's GUI
    polyscope::show();

    return true;
}



#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <queue>
#include <limits>
#include <cmath>

// Function to compute path and distances
void ComputePathAndDistances(
    const Eigen::MatrixXd& border_connected_1,
    int pointIndex,
    int leftIndex,
    int rightIndex,
    double& distanceLeftToPoint,
    double& distanceLeftToRight,
    std::vector<int>& pathLeftToRight
) {
    // Check if indices are valid
    if (leftIndex == -1 || rightIndex == -1) {
        std::cerr << "Error: Left or Right Index not found for the given point.\n";
        distanceLeftToPoint = -1;
        distanceLeftToRight = -1;
        return;
    }

    // Compute the path from leftIndex to rightIndex
    pathLeftToRight.clear();
    for (int i = leftIndex; i <= rightIndex; ++i) {
        pathLeftToRight.push_back(i);
    }

    // Calculate distances
    distanceLeftToRight = 0.0;
    distanceLeftToPoint = 0.0;

    for (size_t i = 0; i < pathLeftToRight.size() - 1; ++i) {
        int currentIndex = pathLeftToRight[i];
        int nextIndex = pathLeftToRight[i + 1];

        double segmentDistance = (border_connected_1.row(currentIndex) - border_connected_1.row(nextIndex)).norm();
        distanceLeftToRight += segmentDistance;

        if (currentIndex == pointIndex || nextIndex == pointIndex) {
            distanceLeftToPoint = distanceLeftToRight;
        }
    }
}



void splitEdgeAndMaintainFaces(
    Eigen::MatrixXd& V2,                // Vertex positions (Nx3)
    Eigen::MatrixXi& F2,                // Faces (Mx3)
    Eigen::MatrixXi& F2_original,       // Original faces (Mx3)
    Eigen::RowVectorXd leftVertex,
    Eigen::RowVectorXd newVertex,
    Eigen::RowVectorXd rightVertex,
    int& info
) {
    // Find indices of leftVertex, newVertex, and rightVertex in V2
    int leftIndex = -1, newIndex = -1, rightIndex = -1;

    for (int i = 0; i < V2.rows(); ++i) {
        if (V2.row(i).isApprox(leftVertex)) {
            leftIndex = i;
        }
        if (V2.row(i).isApprox(newVertex)) {
            newIndex = i;
        }
        if (V2.row(i).isApprox(rightVertex)) {
            rightIndex = i;
        }
    }

    // Ensure all vertices were found
    if (leftIndex == -1 || newIndex == -1 || rightIndex == -1) {
        std::cerr << "Error: One or more vertices not found in V2!" << std::endl;
        return;
    }

    // Debug: Print found indices
    std::cout << "Indices - Left: " << leftIndex << ", New: " << newIndex << ", Right: " << rightIndex << std::endl;

    // Determine the triangles containing the edge (leftIndex, rightIndex)
    std::vector<int> trianglesContainingEdge;
    for (int i = 0; i < F2_original.rows(); ++i) {
        bool containsLeft = (F2_original(i, 0) == leftIndex || F2_original(i, 1) == leftIndex || F2_original(i, 2) == leftIndex);
        bool containsRight = (F2_original(i, 0) == rightIndex || F2_original(i, 1) == rightIndex || F2_original(i, 2) == rightIndex);
        if (containsLeft && containsRight) {
            trianglesContainingEdge.push_back(i);
        }
    }

    // Debug: Print triangles containing the edge
    std::cout << "Triangles containing the edge (" << leftIndex << ", " << rightIndex << "): ";
    for (int t : trianglesContainingEdge) {
        std::cout << t << " ";
    }
    std::cout << std::endl;

    // Function to determine if a triangle is CCW
    auto isTriangleCCW = [&](int A, int B, int C) -> bool {
        Eigen::RowVector3d v0 = V2.row(B) - V2.row(A);
        Eigen::RowVector3d v1 = V2.row(C) - V2.row(A);
        Eigen::RowVector3d normal = v0.cross(v1); // Compute normal
        return normal.z() > 0; // Assumes winding is in the XY plane
    };

    for (int t : trianglesContainingEdge) {
        // Extract the vertices of the triangle
        int A = F2_original(t, 0);
        int B = F2_original(t, 1);
        int C = F2_original(t, 2);

        // Determine indices in the triangle
        int indexLeft = -1, indexRight = -1, indexOther = -1;

        if (A == leftIndex) indexLeft = 0;
        if (B == leftIndex) indexLeft = 1;
        if (C == leftIndex) indexLeft = 2;

        if (A == rightIndex) indexRight = 0;
        if (B == rightIndex) indexRight = 1;
        if (C == rightIndex) indexRight = 2;

        if (A != leftIndex && A != rightIndex) indexOther = 0;
        if (B != leftIndex && B != rightIndex) indexOther = 1;
        if (C != leftIndex && C != rightIndex) indexOther = 2;

        // Debug: Print indices
        std::cout << "Triangle " << t << " - Indices: Left = " << indexLeft
            << ", Right = " << indexRight << ", Other = " << indexOther << std::endl;

        // Ensure indices are valid
        if (indexLeft == -1 || indexRight == -1 || indexOther == -1) {
            std::cerr << "Error: Indices could not be resolved!" << std::endl;
            continue;
        }

        // Create a new triangle ordering to replace the current one
        Eigen::Vector3i updatedTriangle;
        updatedTriangle[indexLeft] = newIndex;
        updatedTriangle[indexRight] = rightIndex;
        updatedTriangle[indexOther] = F2_original(t, indexOther);

        // Replace the current triangle
        //F2.row(t) = updatedTriangle;
        std::cout << "Updated Triangle " << t << ": (" << updatedTriangle[0] << ", "
            << updatedTriangle[1] << ", " << updatedTriangle[2] << ")" << std::endl;

        // Add a new triangle with the same order
        Eigen::Vector3i newTriangle;
        newTriangle[indexLeft] = leftIndex;
        newTriangle[indexRight] = newIndex;
        newTriangle[indexOther] = F2_original(t, indexOther);

        F2.conservativeResize(F2.rows() + 2, Eigen::NoChange);
        F2.row(F2.rows() - 1) = newTriangle;
        F2.row(F2.rows() - 2) = updatedTriangle;

        std::cout << "Added Triangle: (" << newTriangle[0] << ", " << newTriangle[1]
            << ", " << newTriangle[2] << ")" << std::endl;

        //break;
    }
    if (trianglesContainingEdge.size() > 0) {
        info = 1;
    }
    else {
        info = 2;
    }
}









void CompleteBorderCorrespondence(
    Eigen::MatrixXd& V1, Eigen::MatrixXi& F1,
    Eigen::MatrixXd& V2, Eigen::MatrixXi& F2,
    Eigen::MatrixXd& border_1, Eigen::MatrixXd& border_connected_1,
    Eigen::MatrixXd& border_2, Eigen::MatrixXd& border_connected_2
) {
    std::cout << "1,    " << border_1.rows() << " ||| " << border_connected_1.rows() << std::endl;
    std::cout << "2,    " << border_2.rows() << " ||| " << border_connected_2.rows() << std::endl;


    Eigen::MatrixXi F2_original = F2;
    Eigen::MatrixXd border_2_original = border_2;
    

    // Step 1: Find new points in `border_connected_1` that are not in `border_1`
    std::vector<int> newPointsIndices;
    for (int i = 0; i < border_connected_1.rows(); ++i) {
        bool found = false;
        for (int j = 0; j < border_1.rows(); ++j) {
            if (border_connected_1.row(i).isApprox(border_1.row(j), 1e-6)) {
                found = true;
                break;
            }
        }
        if (!found) {
            newPointsIndices.push_back(i);
        }
    }
    Eigen::VectorXi newPointsIndicesEigen = Eigen::Map<Eigen::VectorXi>(newPointsIndices.data(), newPointsIndices.size());
    std::cout << "New indices: " << newPointsIndices.size() << std::endl;

    std::cout << "newPointsIndices: ";
    for (const auto& idx : newPointsIndices) {
        std::cout << idx << " ";
    }
    std::cout << std::endl;

    //polyscope::init();
    //polyscope::registerSurfaceMesh("------ old V2 ", V2, F2);
    //polyscope::registerPointCloud("------- old border V2", border_2);
    

    // Step 2: For each new point, find its leftmost and rightmost neighbors in `border_1`
    for (int idx : newPointsIndices) {
        Eigen::RowVectorXd newPoint = border_connected_1.row(idx);


        std::cout << "----------Part `: Get Locations and path in Mesh1" << std::endl;


        // Find the closest point to the left in `border_1`
        int leftIndex = -1;
        int leftIndexInBorder = -1; // Index in border_1
        for (int i = (idx - 1 + border_connected_1.rows()) % border_connected_1.rows(); i != idx; i = (i - 1 + border_connected_1.rows()) % border_connected_1.rows()) {
            for (int j = 0; j < border_1.rows(); ++j) {
                if (border_connected_1.row(i).isApprox(border_1.row(j), 1e-6)) {
                    leftIndex = i;
                    leftIndexInBorder = j;
                    break;
                }
            }
            if (leftIndex != -1) break; // Exit if left point is found
        }

        // Find the closest point to the right in `border_1`
        int rightIndex = -1;
        int rightIndexInBorder = -1; // Index in border_1
        for (int i = (idx + 1) % border_connected_1.rows(); i != idx; i = (i + 1) % border_connected_1.rows()) {
            for (int j = 0; j < border_1.rows(); ++j) {
                if (border_connected_1.row(i).isApprox(border_1.row(j), 1e-6)) {
                    rightIndex = i;
                    rightIndexInBorder = j;
                    break;
                }
            }
            if (rightIndex != -1) break; // Exit if right point is found
        }


        // Debugging print
        std::cout << "New Point Index: " << idx
            << ", Left Point Index: " << leftIndex
            << ", Right Point Index: " << rightIndex << std::endl;

        if (leftIndex != -1 && rightIndex != -1) {
            std::cout << "Left Point: " << border_connected_1.row(leftIndex)
                << ", Right Point: " << border_connected_1.row(rightIndex) << std::endl;

            

            std::cout << "Left Point Border2: " << border_2_original.row(leftIndexInBorder)
                << ", Right Point Border2: " << border_2_original.row(rightIndexInBorder) << std::endl;

            // Output the corresponding indices in border_1
            std::cout << "Left Point Index in border_1: " << leftIndexInBorder
                << ", Right Point Index in border_1: " << rightIndexInBorder << std::endl;

            // Step 3: Go over the path from left to right and calculate distances
            std::vector<int> pathIndices;
            double distanceLeftToRight = 0.0;
            Eigen::RowVectorXd previousPoint = border_connected_1.row(leftIndex);
            pathIndices.push_back(leftIndex);

            for (int i = leftIndex; i != (rightIndex + 1) % border_connected_1.rows(); i = (i + 1) % border_connected_1.rows()) {
                Eigen::RowVectorXd currentPoint = border_connected_1.row(i);
                pathIndices.push_back(i);

                // Calculate Euclidean distance between consecutive points
                distanceLeftToRight += (currentPoint - previousPoint).norm();
                previousPoint = currentPoint;
            }


            // Calculate distance from left to the new point
            double distanceLeftToPoint = 0.0;
            previousPoint = border_connected_1.row(leftIndex);

            for (int i = leftIndex + 1; i <= idx; ++i) {
                Eigen::RowVectorXd currentPoint = border_connected_1.row(i);
                distanceLeftToPoint += (currentPoint - previousPoint).norm();
                previousPoint = currentPoint;
            }


            // Print results
            std::cout << "Path from Left to Right: ";
            for (int idx : pathIndices) {
                std::cout << idx << " ";
            }
            std::cout << std::endl;


            

            double percentage_distance = std::round(distanceLeftToPoint / distanceLeftToRight * 10e3) / 10e3;

            std::cout << "Distance percentage Left to right : " << percentage_distance << std::endl;
            //std::cout << "Distance from Left to Right: " << distanceLeftToRight << std::endl;
            //std::cout << "Distance from Left to New Point: " << distanceLeftToPoint << std::endl;

            std::cout << "----------Part 2: go over path in second mesh" << std::endl;




            // Step 4: Insert new vertex in V2 between Left Point and Right Point
            // Find the corresponding indices of the left and right points in border_connected_2
            int leftIndexInConnected = -1;
            int rightIndexInConnected = -1;



            for (int i = 0; i < border_connected_2.rows(); ++i) {
                //std::cout << i << std::endl;
                if (border_connected_2.row(i).isApprox(border_2_original.row(leftIndexInBorder), 1e-6)) {
                    leftIndexInConnected = i;
                }
                if (border_connected_2.row(i).isApprox(border_2_original.row(rightIndexInBorder), 1e-6)) {
                    rightIndexInConnected = i;
                }


    

               
                if (leftIndexInConnected != -1 && rightIndexInConnected != -1 ) break;

            }

            std::cout << "in border_2_connected, left, rightindex " << leftIndexInConnected << " ||| " << rightIndexInConnected << std::endl;
             
            // now you got the original index in the border ones. 
            // use these to then find the index in the updated border_2


            if (leftIndexInConnected == -1 || rightIndexInConnected == -1 ) {
                //std::cerr << "Error: Could not find left or right index in border_connected_2!" << std::endl;
                std::cout << "------------------------------------ vertex doesn't exist (THIS SHOULDNt happen)" << std::endl;

                continue;
            }

            //edge case lefindexInconnected other way around!?!
            // so swap them then. 
            // Swap the values -- might not be bugfree, but highly likely no issues
            if (leftIndexInConnected > rightIndexInConnected && rightIndexInConnected != 0) {
                std::swap(leftIndexInConnected, rightIndexInConnected);
            }


            // Compute the total length of edges between leftIndexInConnected and rightIndexInConnected
            double totalEdgeLength = 0.0;
            int currentIndex = leftIndexInConnected;

            while (true) {
                int nextIndex = (currentIndex + 1) % border_connected_2.rows();

                // Break condition: stop after the edge that ends at rightIndexInConnected
                if (currentIndex == rightIndexInConnected) break;


                Eigen::RowVectorXd currentPoint = border_connected_2.row(currentIndex);
                Eigen::RowVectorXd nextPoint = border_connected_2.row(nextIndex);
                //std::cout << "index first: " << currentIndex << " point: " << currentPoint << std::endl;
                //std::cout << "index second: " << nextIndex << " point: " << nextPoint << std::endl;

                totalEdgeLength += (nextPoint - currentPoint).norm();
                currentIndex = nextIndex; // Move to the next index
            }

            //std::cout << "Total Edge Length: " << totalEdgeLength << std::endl;

            // Calculate the target distance based on percentage_distance
            double targetDistance = percentage_distance * totalEdgeLength;
            //std::cout << "Target Distance (Percentage): " << targetDistance << std::endl;

            // Walk along the edges to reach the target distance
            double accumulatedDistance = 0.0;
            currentIndex = leftIndexInConnected; // Start from the left index

            while (true) {
                int nextIndex = (currentIndex + 1) % border_connected_2.rows();

                // Stop when the range ends
                if (currentIndex == rightIndexInConnected) break;

                Eigen::RowVectorXd currentPoint = border_connected_2.row(currentIndex);
                Eigen::RowVectorXd nextPoint = border_connected_2.row(nextIndex);
                double edgeLength = (nextPoint - currentPoint).norm();

                accumulatedDistance += edgeLength;

                //std::cout << "Walking edge: " << currentIndex << " -> " << nextIndex << std::endl;
                //std::cout << "  Current Point: " << currentPoint.transpose() << std::endl;
                //std::cout << "  Previous Point: " << border_connected_2.row(currentIndex).transpose() << std::endl;
                //std::cout << "  Edge Length: " << edgeLength << std::endl;
                //std::cout << "  Accumulated Distance: " << accumulatedDistance << std::endl;

                if (accumulatedDistance >= targetDistance) {
                    //std::cout << "  Found insertion edge at index: " << currentIndex << std::endl;
                    break;
                }

                currentIndex = nextIndex; // Move to the next index
            }

            // Calculate the exact position of the new vertex along the identified edge
            Eigen::RowVectorXd edgeStart = border_connected_2.row(currentIndex);
            Eigen::RowVectorXd edgeEnd = border_connected_2.row((currentIndex + 1) % border_connected_2.rows());
            double edgeLength = (edgeEnd - edgeStart).norm();
            double edgePercentage = (targetDistance - (accumulatedDistance - edgeLength)) / edgeLength;

            Eigen::RowVectorXd newVertex = edgeStart + (edgeEnd - edgeStart) * edgePercentage;

            //std::cout << "Final Calculations:" << std::endl;
            //std::cout << "  Edge Start: " << edgeStart.transpose() << std::endl;
            //std::cout << "  Edge End: " << edgeEnd.transpose() << std::endl;
            //std::cout << "  Edge Length: " << edgeLength << std::endl;
            //std::cout << "  Edge Percentage: " << edgePercentage << std::endl;
            //std::cout << "  New Vertex: " << newVertex.transpose() << std::endl;


            Eigen::MatrixXd V2_backup = V2;

            // Insert the new vertex into V2
            V2.conservativeResize(V2.rows() + 1, Eigen::NoChange);
            V2.row(V2.rows() - 1) = newVertex;





            // now given edgeStart, edgeEnd and newVertex on that line between those points.
            // split the faces ath the newVertex location. Every edge is used by 2 triangles. Split them into 4
            // so triangle 1 : A, B, C, tirangle 2, A, B, D.  Where A=edgestart, B=edgeEnd (or B=edgeStart, A=edgeEnd). X = newVertex Split this into
            // T1: A,X,C   
            // T2: X,B,C
            // T3: A, X, D
            // T4: X, B, D

            //input for function: V2, F2, (A,B,X)

            // make sure the references of the faces maintain ok..
            int infoX = 0;
            splitEdgeAndMaintainFaces(
                V2,                // Vertex positions (Nx3)
                F2,                // Faces (Mx3)
                F2_original, 
                edgeStart,
                newVertex,
                edgeEnd,
                infoX
            );
            if (infoX == 1 || infoX == 2) {

                // Update border_2 by inserting the new vertex at the appropriate position
                Eigen::MatrixXd updatedBorder2(border_2.rows() + 1, border_2.cols());
                for (int i = 0; i < idx; ++i) {
                    updatedBorder2.row(i) = border_2.row(i);
                }
                updatedBorder2.row(idx) = newVertex;
                for (int i = idx; i < border_2.rows(); ++i) {
                    updatedBorder2.row(i + 1) = border_2.row(i);
                }
                border_2 = updatedBorder2;



                std::cout << "New vertex added at index: " << idx << "Vertex:" << newVertex<<  std::endl;


                //if (idx == 366 || idx == 757 || idx == 1747) {
                //    std::string name;

                //    std::getline(std::cin, name);  // To take a full line input, including spaces
                //}
            }
            if (infoX == 2 ){
                //std::cout << " indx could be faulty   =======================  " << idx <<  std::endl;
                //    std::string x;
                //    std::cin >> x;
                
                //V2 = V2_backup;
            }

        }


        std::cout << "-----------------\n\n" << std::endl;



    }

    std::cout << "1 - updated,    " << border_1.rows() << " ||| " << border_connected_1.rows() << std::endl;
    std::cout << "2 - updated,    " << border_2.rows() << " ||| " << border_connected_2.rows() << std::endl;


    //polyscope::registerSurfaceMesh("------ new V2 ", V2, F2);
    //polyscope::registerPointCloud("------- new border V2", border_2);

    //polyscope::show();

}







