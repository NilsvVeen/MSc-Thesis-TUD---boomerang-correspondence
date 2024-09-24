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
        viewSTLObject(modelPath);

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
        viewSTLObject(modelPath);

        // get outlining vertices of object
        std::vector<Eigen::Vector2d> border_vertices = fitPlaneAndAlignMesh(modelPath, directoryName2);
    }



    // view models and parameterize
    if (ParameterizeObjects) {

        polyscope::removeAllGroups();
        polyscope::removeAllStructures();

        Eigen::MatrixXd V;
        Eigen::MatrixXd V2;

        // Reading point cloud from PLY file
        if (readPointCloudFromFile(directoryName + "/alpha_shape_border.obj", V) 
            && readPointCloudFromFile(directoryName + "/alpha_shape_border.obj", V2)
            ) {
            std::cout << "Point cloud loaded successfully." << std::endl;
            std::cout << "Number of vertices: " << V.rows() << std::endl; // Print the number of vertices
            std::cout << "Number of vertices in 2: " << V2.rows() << std::endl; // Print the number of vertices


            // Sort vertices by proximity
            Eigen::MatrixXd sortedVertices = sortVerticesByProximity(V);
            savePointCloudToFile(directoryName + "/border_vertices_in_order.obj", sortedVertices);

            Eigen::MatrixXd sortedVertices2 = sortVerticesByProximity(V2);
            savePointCloudToFile(directoryName2 + "/border_vertices_in_order.obj", sortedVertices2);

            //showPointCloudInParts(sortedVertices);

            //std::vector<double> unit_parameters = computeUnitParametrization(sortedVertices);

            //writeParamsToFile("output_files/unit_parameters.txt", unit_parameters);

            //showSelection(sortedVertices);

            showSideBySideSelectionWithVertexSelection(sortedVertices, sortedVertices2);
        }
    }

    return 0;
}



