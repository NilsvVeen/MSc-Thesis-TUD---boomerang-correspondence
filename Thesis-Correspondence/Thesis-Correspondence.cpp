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
    std::string directoryName = "output_files";
    createDirectory(directoryName);

    if (ProcessObjects) {
        old_1(directoryName);
    }

   
    if (ProcessObjects) {

        // Print files in the global models directory
        printFilesInDirectory(GLOBAL_MODELS_DIRECTORY);

        std::string modelPath = GLOBAL_MODELS_DIRECTORY + "/Boomerang_02.stl";

        // Call the function to view the STL object
        viewSTLObject(modelPath);

        // Call the function to compute and save the normal direction
        std::vector<Eigen::Vector2d> border_vertices = fitPlaneAndAlignMesh(modelPath, directoryName);
    }


    if (ParameterizeObjects) {

        Eigen::MatrixXd V;

        // Reading point cloud from PLY file
        if (readPointCloudFromFile(directoryName + "/alpha_shape_border.obj", V)) {
            std::cout << "Point cloud loaded successfully." << std::endl;
            std::cout << "Number of vertices: " << V.rows() << std::endl; // Print the number of vertices


            // Sort vertices by proximity
            Eigen::MatrixXd sortedVertices = sortVerticesByProximity(V);
            savePointCloudToFile("output_files/border_vertices_in_order.obj", sortedVertices);

            //showPointCloudInParts(sortedVertices);

            //std::vector<double> unit_parameters = computeUnitParametrization(sortedVertices);

            //writeParamsToFile("output_files/unit_parameters.txt", unit_parameters);

            //showSelection(sortedVertices);

            showSideBySideSelectionWithVertexSelection(sortedVertices, sortedVertices);
        }
    }

    return 0;
}



