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





