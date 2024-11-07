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
#include "splitMesh.h"

const std::string GLOBAL_MODELS_DIRECTORY = MODELS_DIRECTORY;

static const bool ProcessObjects = false;
static const bool ParameterizeObjects = false;
static const bool ReadCalculateSortedVertices = false; // only enable if already calculated
static const bool showOriginalRotatedMesh = false;
static const bool correspondences2dto3d = false;
static const bool parameterizeSurfaceBool = true;


int main()
{

    //// Create a directory for output files
    std::string directoryName = "output_files_3";
    createDirectory(directoryName);
    std::string directoryName2 = "output_files_4";
    createDirectory(directoryName2);


    // get 2d outline
    if (ProcessObjects) {

        polyscope::removeAllGroups();
        polyscope::removeAllStructures();

        // Print files in the global models directory
        printFilesInDirectory(GLOBAL_MODELS_DIRECTORY);

        std::string modelPath = GLOBAL_MODELS_DIRECTORY + "/Boomerang_09.stl";

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

        std::string modelPath = GLOBAL_MODELS_DIRECTORY + "/Boomerang_12.stl";

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

        if (true) {

            std::cout << "split Meshes" << std::endl;

            const std::string splitmesh = "SplitMeshes";
            createDirectory(splitmesh);
            clearDirectory(splitmesh);


            Eigen::MatrixXd MeshA_V;
            Eigen::MatrixXi MeshA_F;
            Eigen::MatrixXd MeshB_V;
            Eigen::MatrixXi MeshB_F;

            
            std::cout << "Step 1: remove lifted curve vertices and faces" << std::endl;
            RemoveVerticesAndFaces(Mesh1_V, Mesh1_F, V3, MeshA_V, MeshA_F);
            countConnectedComponents(MeshA_F);

            Eigen::MatrixXd border_V = getBorderVerticesMatrix(MeshA_V, MeshA_F);

            //showMeshAndPointCloud(MeshA_V, MeshA_F, border_V);


            Eigen::MatrixXd removed_V;
            Eigen::MatrixXi removed_F;
            Eigen::MatrixXd border_V_new;


            std::cout << "Step 2: remove leftover Bridges" << std::endl;
            removeVerticesWithTwoFacesAndBorderEdges(MeshA_V, MeshA_F, border_V, MeshB_V, MeshB_F, removed_V, removed_F, border_V_new);
            countConnectedComponents(MeshB_F);

            //showMeshAndPointCloud(MeshB_V, MeshB_F, border_V_new);


            Eigen::MatrixXd Mesh_V_Split1;
            Eigen::MatrixXi Mesh_F_Split1;
            Eigen::MatrixXd Mesh_V_Split2;
            Eigen::MatrixXi Mesh_F_Split2;
            //split into 2 meshes
            splitMeshIn2(MeshB_V,
                MeshB_F, Mesh_V_Split1, Mesh_F_Split1, Mesh_V_Split2, Mesh_F_Split2);


            std::cout << "revmoedV" << std::endl;
            std::cout << removed_V << std::endl;
            std::cout << "removedF" << std::endl;

            std::cout << removed_F << std::endl;

            std::cout << "----------------- Add back faces to 1" << std::endl;
            //AddBackFace(Mesh_V_Split1, Mesh_F_Split1, removed_V, removed_F);
            FindMatchingEdges(MeshA_V, MeshA_F, removed_V, removed_F, 1e-4, Mesh_V_Split1, Mesh_F_Split1);
            FindMatchingEdges(MeshA_V, MeshA_F, removed_V, removed_F, 1e-4, Mesh_V_Split2, Mesh_F_Split2);

            polyscope::init();
            polyscope::registerSurfaceMesh("Mesh INput", Mesh1_V, Mesh1_F);
            polyscope::registerSurfaceMesh("Mesh AAA", MeshA_V, MeshA_F);
            polyscope::registerSurfaceMesh("Mesh BBB", MeshB_V, MeshB_F);
            polyscope::registerSurfaceMesh("Mesh P1", Mesh_V_Split1, Mesh_F_Split1);
            polyscope::registerSurfaceMesh("Mesh P2", Mesh_V_Split2, Mesh_F_Split2);
            polyscope::registerPointCloud("INput lifted curve", V3);
            polyscope::registerPointCloud("Removed A to B", removed_V);
            polyscope::registerPointCloud("Border A to B", border_V_new);
            polyscope::show();


            //saveMeshToFile(splitmesh + "/A.obj", Mesh_V_Split1, Mesh_F_Split1);
            //saveMeshToFile(splitmesh + "/B.obj", Mesh_V_Split2, Mesh_F_Split2);


            //// parameterize whole surface
            //if (true) {
            //    Eigen::MatrixXd UV_split1;
            //    Eigen::MatrixXd V_border_split1 = getBorderVerticesMatrix(Mesh_V_Split1, Mesh_F_Split1);
            //    if (!paramsurface5(Mesh_V_Split1, Mesh_F_Split1, UV_split1, V_border_split1, false)) {
            //        std::cerr << "Surface parameterization failed.\n";
            //        return EXIT_FAILURE;
            //    }
            //}
            //// parameterize whole surface
            //if (true) {
            //    Eigen::MatrixXd UV_split2;
            //    Eigen::MatrixXd V_border_split2 = getBorderVerticesMatrix(Mesh_V_Split2, Mesh_F_Split2);
            //    if (!paramsurface5(Mesh_V_Split2, Mesh_F_Split2, UV_split2, V_border_split2, false)) {
            //        std::cerr << "Surface parameterization failed.\n";
            //        return EXIT_FAILURE;
            //    }
            //}

        }


        // parameterize whole surface
        if (false) {
            Eigen::MatrixXd UV1;
            if (!paramsurface5(Mesh1_V, Mesh1_F, UV1, V3, true)) {
                std::cerr << "Surface parameterization failed.\n";
                return EXIT_FAILURE;
            }
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



