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
#include "2dCorrespondenceTo3d.h"

#include "parameterizeSurface.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif
#include <regex>

const std::string GLOBAL_MODELS_DIRECTORY = MODELS_DIRECTORY;

static const bool CircleElipFlag = false;
static const bool ProcessObjects = false;
static const bool ParameterizeObjects = false;
static const bool ReadCalculateSortedVertices = false; // only enable if already calculated
static const bool showOriginalRotatedMesh = false;
static const bool correspondences2dto3d = false;
static const bool parameterizeSurfaceBool = true;


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

    {

    //if (false) {
    //    // Load the first mesh
    //    Eigen::MatrixXd Mesh1_V;
    //    Eigen::MatrixXi Mesh1_F;
    //    readMeshFromFile("temp_test/boom2.obj", Mesh1_V, Mesh1_F);

    //    // Save the mesh to another file
    //    saveMeshToFile("temp_test/boom2_saved.obj", Mesh1_V, Mesh1_F);

    //    // Load the second mesh from the saved file
    //    Eigen::MatrixXd Mesh2_V;
    //    Eigen::MatrixXi Mesh2_F;
    //    readMeshFromFile("temp_test/boom2_saved.obj", Mesh2_V, Mesh2_F);

    //    // Print information about the meshes
    //    std::cout << "Mesh 1:\n";
    //    std::cout << "Vertices: " << Mesh1_V.rows() << " x " << Mesh1_V.cols() << "\n";
    //    std::cout << "Faces: " << Mesh1_F.rows() << " x " << Mesh1_F.cols() << "\n\n";

    //    std::cout << "Mesh 2:\n";
    //    std::cout << "Vertices: " << Mesh2_V.rows() << " x " << Mesh2_V.cols() << "\n";
    //    std::cout << "Faces: " << Mesh2_F.rows() << " x " << Mesh2_F.cols() << "\n";

    //    // Check if Mesh1 and Mesh2 are identical and print differences
    //    bool identical = true;

    //    double tolerance = 1e-3; // Custom tolerance for checking approximate equality

    //    // Check vertices
    //    if (!Mesh1_V.isApprox(Mesh2_V, tolerance)) {
    //        identical = false;
    //        std::cout << "\nDifferences found in vertices (greater than " << tolerance << "):\n";
    //        for (int i = 0; i < Mesh1_V.rows(); ++i) {
    //            for (int j = 0; j < Mesh1_V.cols(); ++j) {
    //                if (std::abs(Mesh1_V(i, j) - Mesh2_V(i, j)) > tolerance) {
    //                    std::cout << "Vertex difference at (" << i << ", " << j << "): "
    //                        << std::setprecision(10)
    //                        << "Mesh1_V = " << Mesh1_V(i, j) << ", Mesh2_V = " << Mesh2_V(i, j) << "\n";
    //                }
    //            }
    //        }
    //    }


    //    // Check faces
    //    if (!Mesh1_F.isApprox(Mesh2_F)) {
    //        identical = false;
    //        std::cout << "\nDifferences found in faces:\n";
    //        for (int i = 0; i < Mesh1_F.rows(); ++i) {
    //            for (int j = 0; j < Mesh1_F.cols(); ++j) {
    //                if (Mesh1_F(i, j) != Mesh2_F(i, j)) {
    //                    std::cout << "Face difference at (" << i << ", " << j << "): "
    //                        << "Mesh1_F = " << Mesh1_F(i, j) << ", Mesh2_F = " << Mesh2_F(i, j) << "\n";
    //                }
    //            }
    //        }
    //    }

    //    if (identical) {
    //        std::cout << "\nThe meshes are identical.\n";
    //    }
    //    else {
    //        std::cout << "\nThe meshes differ.\n";
    //    }


    //    // Load the second mesh from the saved file
    //    Eigen::MatrixXd Mesh3_V;
    //    Eigen::MatrixXi Mesh3_F;
    //    readMeshFromFile("output_files_1/rotated_mesh.obj", Mesh3_V, Mesh3_F);

    //    std::cout << "Mesh 3:\n";
    //    std::cout << "Vertices: " << Mesh3_V.rows() << " x " << Mesh3_V.cols() << "\n";
    //    std::cout << "Faces: " << Mesh3_F.rows() << " x " << Mesh3_F.cols() << "\n";



    //}

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


    // Create default directory if none is provided
    const std::string DEFAULT_CORRESPONDENCES_meshes_FOLDER = "models_for_correspondences";

    if (showOriginalRotatedMesh) {
        polyscope::init();
        polyscope::removeAllGroups();
        polyscope::removeAllStructures();

        Eigen::MatrixXd Mesh1_V;
        Eigen::MatrixXi Mesh1_F;
        readMeshFromFile(directoryName + "/rotated_mesh.obj", Mesh1_V, Mesh1_F);

        Eigen::MatrixXd Mesh2_V;
        Eigen::MatrixXi Mesh2_F;
        readMeshFromFile(directoryName2 + "/rotated_mesh.obj", Mesh2_V, Mesh2_F);

        Eigen::MatrixXd V1_border;
        Eigen::MatrixXd V2_border;
        readPointCloudFromFile(directoryName + "/alpha_shape_border.obj", V1_border);
        readPointCloudFromFile(directoryName2 + "/alpha_shape_border.obj", V2_border);

        // Display the two meshes side by side
        showSideBySideMeshes(Mesh1_V, Mesh1_F, Mesh2_V, Mesh2_F, V1_border, V2_border, DEFAULT_CORRESPONDENCES_meshes_FOLDER);
    }



    const std::string DEFAULT_CORRESPONDENCES_FOLDER = "Correspondences";
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

            showSideBySideSelectionWithVertexSelection(sortedVertices, sortedVertices2, DEFAULT_CORRESPONDENCES_FOLDER);



            //showSideBySideSelectionWithVertexSelection(sampled.first, sampled.second);





        }
    }

    const std::string DEFAULT_2dto3d_FOLDER = "2d_Curve_in_3d";
    if (correspondences2dto3d) {
        const std::string meshesFolder = DEFAULT_CORRESPONDENCES_meshes_FOLDER;
        const std::string pointCloudsFolder = DEFAULT_CORRESPONDENCES_FOLDER;
        const std::string Curve2dTo3dFolder = DEFAULT_2dto3d_FOLDER;

        // Declare variables to hold the data
        Eigen::MatrixXd Mesh1_V, Mesh2_V;
        Eigen::MatrixXi Mesh1_F, Mesh2_F;
        std::vector<Eigen::MatrixXd> V1_pointclouds;
        std::vector<Eigen::MatrixXd> V2_pointclouds;

        // Call the function to read the data directly into the variables
        readMeshesAndPointClouds(meshesFolder, pointCloudsFolder,
            Mesh1_V, Mesh1_F,
            Mesh2_V, Mesh2_F,
            V1_pointclouds, V2_pointclouds);


        std::cout << "match V1" << std::endl;
        // Process correspondences for V1 (Exact match x, y, use z from Mesh1)
        for (Eigen::MatrixXd& V1 : V1_pointclouds) {
            findExactCorrespondences(Mesh1_V, V1);
        }

        std::cout << "match V2" << std::endl;

        // Process correspondences for V2 (Closest x, y match, use x, y, z from Mesh2)
        for (Eigen::MatrixXd& V2 : V2_pointclouds) {
            findClosestCorrespondences(Mesh2_V, V2);
        }


        // Check if meshes were successfully read
        if (Mesh1_V.rows() == 0 || Mesh2_V.rows() == 0) {
            std::cerr << "Error: Meshes not successfully loaded." << std::endl;
            return -1;
        }

        writeOutputsToFolder(Curve2dTo3dFolder, Mesh1_V, Mesh1_F, Mesh2_V, Mesh2_F, V1_pointclouds, V2_pointclouds);
 
        // Show them in Polyscope with the common color
        showInPolyscope(Mesh1_V, Mesh1_F, Mesh2_V, Mesh2_F, V1_pointclouds, V2_pointclouds);

    }


    if (parameterizeSurfaceBool) {

        std::cout << "surface Parameterization:" << std::endl;

        const std::string surfaceParam = "surfaceParameterize";
        createDirectory(surfaceParam);
        clearDirectory(surfaceParam);

        Eigen::MatrixXd Mesh1_V;
        Eigen::MatrixXi Mesh1_F;
        readMeshFromFile( "2d_Curve_in_3d/Mesh1.obj", Mesh1_V, Mesh1_F);


        std::string V1_regex = "V1_pointcloud_(\\d+)\\.txt";
        std::string V2_regex = "V1_pointcloud_(\\d+)\\.txt";

        Eigen::MatrixXd V3 = readAndConcatenatePointClouds(DEFAULT_2dto3d_FOLDER, V1_regex);


        Eigen::MatrixXd UV1;
        if (!paramsurface5(Mesh1_V, Mesh1_F, UV1, V3)) {
            std::cerr << "Surface parameterization failed.\n";
            return EXIT_FAILURE;
        }



        //saveMeshToFile(DEFAULT_CORRESPONDENCES_meshes_FOLDER + "UV.obj", UV1, Mesh1_F);



        //Eigen::MatrixXd Mesh2_V;
        //Eigen::MatrixXi Mesh2_F;
        //readMeshFromFile(DEFAULT_CORRESPONDENCES_meshes_FOLDER + "/RigthMesh.obj", Mesh2_V, Mesh2_F);

        //Eigen::MatrixXd V1_border;
        //Eigen::MatrixXd V2_border;
        //readPointCloudFromFile(directoryName + "/alpha_shape_border.obj", V1_border);
        //readPointCloudFromFile(directoryName2 + "/alpha_shape_border.obj", V2_border);

        //// Display the two meshes side by side
        //showSideBySideMeshes(Mesh1_V, Mesh1_F, Mesh2_V, Mesh2_F, V1_border, V2_border, DEFAULT_CORRESPONDENCES_meshes_FOLDER);
    }



    return 0;
}



