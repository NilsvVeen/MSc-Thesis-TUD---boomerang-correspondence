// Thesis-Correspondence.cpp : Defines the entry point for the application.
//

#include "Thesis-Correspondence.h"

// main.cpp
#include "geometry.h"
#include "file_utils.h"
#include "sfml_utils.h"
#include "stl_utils.h"
#include "border_vertex_extraction.h"
#include "stl_utils.h"
#include "parameterize.h"
#include <polyscope/polyscope.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

const std::string GLOBAL_MODELS_DIRECTORY = MODELS_DIRECTORY;

static const bool CircleElipFlag = false;
static const bool ProcessObjects = false;
static const bool ParameterizeObjects = true;
static const bool ReadCalculateSortedVertices = true;
static const bool showOriginalRotatedMesh = true;

int main()
{

    //// Create a directory for output files
    std::string directoryName = "output_files_1";
    createDirectory(directoryName);
    std::string directoryName2 = "output_files_2";
    createDirectory(directoryName2);

    // unused
    if (CircleElipFlag) {
        old_1(directoryName);
    }

   
    // get 2d outline
    if (ProcessObjects) {

        polyscope::removeAllGroups();
        polyscope::removeAllStructures();

        // Print files in the global models directory
        printFilesInDirectory(GLOBAL_MODELS_DIRECTORY);

        std::string modelPath = GLOBAL_MODELS_DIRECTORY + "/Boomerang_02.stl";

        // Call the function to view the STL object
        //viewSTLObject(modelPath);

        // get outlining vertices of object
        std::vector<Eigen::Vector2d> border_vertices = fitPlaneAndAlignMesh(modelPath, directoryName);
    }

    // get 2d outline
    if (ProcessObjects) {

        polyscope::removeAllGroups();
        polyscope::removeAllStructures();

        // Print files in the global models directory
        printFilesInDirectory(GLOBAL_MODELS_DIRECTORY);

        std::string modelPath = GLOBAL_MODELS_DIRECTORY + "/Boomerang_03.stl";

        // Call the function to view the STL object
        //viewSTLObject(modelPath);

        // get outlining vertices of object
        std::vector<Eigen::Vector2d> border_vertices = fitPlaneAndAlignMesh(modelPath, directoryName2);
    }


    if (showOriginalRotatedMesh) {
        polyscope::init();

        Eigen::MatrixXd Mesh1_V;
        Eigen::MatrixXi Mesh1_F;
        readMeshFromFile(directoryName + "/rotated_mesh.obj", Mesh1_V, Mesh1_F);

        Eigen::MatrixXd Mesh2_V;
        Eigen::MatrixXi Mesh2_F;
        readMeshFromFile(directoryName2 + "/rotated_mesh.obj", Mesh2_V, Mesh2_F);

        // Display the two meshes side by side
        showSideBySideMeshes(Mesh1_V, Mesh1_F, Mesh2_V, Mesh2_F);
    }




    // view models and parameterize
    if (ParameterizeObjects) {

        // clear crap!
        //polyscope::removeAllGroups();
        //polyscope::removeAllStructures();

        Eigen::MatrixXd V;
        Eigen::MatrixXd V2;






        // Reading point cloud from PLY file
        if (readPointCloudFromFile(directoryName + "/alpha_shape_border.obj", V) 
            && readPointCloudFromFile(directoryName2 + "/alpha_shape_border.obj", V2)
            ) {
            std::cout << "Point cloud loaded successfully." << std::endl;
            std::cout << "Number of vertices: " << V.rows() << std::endl; // Print the number of vertices
            std::cout << "Number of vertices in 2: " << V2.rows() << std::endl; // Print the number of vertices


            Eigen::MatrixXd sortedVertices, sortedVertices2;

            if (ReadCalculateSortedVertices) {
                // Read sorted vertices from file
                readPointCloudFromFile(directoryName + "/border_vertices_in_order.obj", sortedVertices);
                readPointCloudFromFile(directoryName2 + "/border_vertices_in_order.obj", sortedVertices2);
            }
            else {
                // Sort vertices by proximity and reverse order
                sortedVertices = reverseOrder(sortVerticesByProximity(V));
                sortedVertices2 = reverseOrder(sortVerticesByProximity(V2));

                // Write sorted vertices to file
                savePointCloudToFile(directoryName + "/border_vertices_in_order.obj", sortedVertices);
                savePointCloudToFile(directoryName2 + "/border_vertices_in_order.obj", sortedVertices2);
            }



            //sample the points so both point clouds have equal amount of points.
            //auto sampled = getEqualizedPointClouds(sortedVertices, sortedVertices);
            //savePointCloudToFile(directoryName + "/border_vertices_in_order_sampled.obj", sampled.first);
            //savePointCloudToFile(directoryName2 + "/border_vertices_in_order_sampled.obj", sampled.second);

            //showPointCloudInParts(sortedVertices);

            std::vector<double> unit_parameters = computeUnitParametrization(sortedVertices);

            writeParamsToFile("output_files/unit_parameters.txt", unit_parameters);

            //showSelection(sortedVertices);

            showSideBySideSelectionWithVertexSelection(sortedVertices, sortedVertices2);



            //showSideBySideSelectionWithVertexSelection(sampled.first, sampled.second);




//// Clockwise circle points (with z = 0)
//            Eigen::MatrixXd V1_test(8, 3);
//            V1_test << 1.0, 0.0, 0.0,    // Rightmost point on the circle
//                0.71, -0.71, 0.0, // Bottom-right
//                0.0, -1.0, 0.0,   // Bottom
//                -0.71, -0.71, 0.0, // Bottom-left
//                -1.0, 0.0, 0.0,    // Leftmost point on the circle
//                -0.71, 0.71, 0.0,  // Top-left
//                0.0, 1.0, 0.0,    // Top
//                0.71, 0.71, 0.0;  // Top-right
//
//            // Anti-clockwise circle points (with z = 0)
//            Eigen::MatrixXd V2_test(8, 3);
//            V2_test << 1.0, 0.0, 0.0,    // Rightmost point on the circle
//                0.71, 0.71, 0.0,  // Top-right
//                0.0, 1.0, 0.0,    // Top
//                -0.71, 0.71, 0.0,  // Top-left
//                -1.0, 0.0, 0.0,    // Leftmost point on the circle
//                -0.71, -0.71, 0.0, // Bottom-left
//                0.0, -1.0, 0.0,   // Bottom
//                0.71, -0.71, 0.0; // Bottom-right
//
//            showSideBySideSelectionWithVertexSelection(V1_test, V1_test);


        }
    }

    return 0;
}



