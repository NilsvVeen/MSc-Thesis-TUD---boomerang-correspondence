// Thesis-Correspondence.cpp : Defines the entry point for the application.
//

#include "Thesis-Correspondence.h"

// main.cpp
#include "geometry.h"
#include "file_utils.h"
#include "stl_utils.h"
#include "border_vertex_extraction.h"
#include "stl_utils.h"
#include "parameterize.h"
#include <polyscope/polyscope.h>
#include "2dCorrespondenceTo3d.h"

#include "parameterizeSurface.h"

//#ifndef PI
//#define PI 3.14159265358979323846
//#endif
#include <regex>
#include "splitMesh.h"
#include "uvMapCorrespondence.h"
#include "evaluateCorrespondence.h"
#include "intermediateShape.h"


const std::string GLOBAL_MODELS_DIRECTORY = MODELS_DIRECTORY;
static const bool ReadCalculateSortedVertices = false; // only enable if already calculated



static const bool ProcessObjects = false;
static const bool ParameterizeObjects = false;


static const bool shiftAll = false;


static const bool correspondences2dto3d = false;

static const bool parameterizeSurfaceBool = false;
static const bool uvMapCorrespondence = false;
static const bool evaluateCorrespondence = false;




static const bool newShapeMake = true;



int main()
{



    //// Create a directory for output files
    std::string directoryName = "output_files_3";
    createDirectory(directoryName);
    std::string directoryName2 = "output_files_4";
    createDirectory(directoryName2);

    std::cout << "Step 1 ---------- Process objects to get boundary vertices" << std::endl;

    // get 2d outline
    if (ProcessObjects) {

        

        polyscope::removeAllGroups();
        polyscope::removeAllStructures();

        // Print files in the global models directory
        printFilesInDirectory(GLOBAL_MODELS_DIRECTORY);

        //std::string modelPath = GLOBAL_MODELS_DIRECTORY + "/Boomerang_09.stl";
        std::string modelPath = GLOBAL_MODELS_DIRECTORY + "/Boomerang_09_decimate01.stl";

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

        //std::string modelPath = GLOBAL_MODELS_DIRECTORY + "/Boomerang_12.stl";
        //std::string modelPath = GLOBAL_MODELS_DIRECTORY + "/Boomerang_12_decimate01.stl";
        //std::string modelPath = GLOBAL_MODELS_DIRECTORY + "/Boomerang_12_decimate01.stl";
        std::string modelPath = GLOBAL_MODELS_DIRECTORY + "/Boomerang_11_decimate01.stl";

        // Call the function to view the STL object
        //viewSTLObject(modelPath);

        Eigen::MatrixXd Mesh1_V;
        Eigen::MatrixXi Mesh1_F;
        Eigen::MatrixXd V1_border;
        readMeshFromFile(directoryName + "/rotated_mesh.obj", Mesh1_V, Mesh1_F);
        readPointCloudFromFile(directoryName + "/alpha_shape_border.obj", V1_border);

        // get outlining vertices of object
        std::vector<Eigen::Vector2d> border_vertices = fitPlaneAndAlignMesh(modelPath, directoryName2);

    }





    std::cout << "Step 3 ---------- sort border vertices and let user select vertices to get correspondences" << std::endl;

    Eigen::Vector3d shiftVector;

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


            showSideBySideSelectionWithVertexSelection(sortedVertices, sortedVertices2, DEFAULT_CORRESPONDENCES_FOLDER, shiftVector);
        }
    }

    std::cout << "official shift vector to use!!!!!!!!!!!!!!!!! : " << shiftVector << std::endl;


    const std::string shiftMeshAndCurve = "shiftMeshAndCurve";
    if (shiftAll) {
        std::cout << "Step 3.5 ---------- Shift Mesh 2 by shift vector" << std::endl;

        //shiftVector << 538, 496 , 349; // normal case
        //shiftVector << 538, 505 , 349; // decimated 01  case
        //shiftVector << 856, 84,  317; // decimated 01  case
        shiftVector << 913 ,84 ,319; // decimated 01  case
        //shiftVector << 213, 0 , 0; // decimated 01  case



        Eigen::RowVector3d shiftVectorRow = shiftVector.transpose();

        createDirectory(shiftMeshAndCurve);
        clearDirectory(shiftMeshAndCurve);

        Eigen::MatrixXd sortedVertices, sortedVertices2_temp;
        readPointCloudFromFile(directoryName + "/border_vertices_in_order.obj", sortedVertices);
        readPointCloudFromFile(directoryName2 + "/border_vertices_in_order.obj", sortedVertices2_temp);
        Eigen::MatrixXd sortedVertices2 = sortedVertices2_temp.rowwise() + shiftVectorRow;

        Eigen::MatrixXd Mesh1_V;
        Eigen::MatrixXi Mesh1_F;
        readMeshFromFile(directoryName + "/rotated_mesh.obj", Mesh1_V, Mesh1_F);

        Eigen::MatrixXd Mesh2_V_temp;
        Eigen::MatrixXi Mesh2_F;
        readMeshFromFile(directoryName2 + "/rotated_mesh.obj", Mesh2_V_temp, Mesh2_F);
        Eigen::MatrixXd Mesh2_V = Mesh2_V_temp.rowwise() + shiftVectorRow;


        saveMeshToFile(shiftMeshAndCurve + "/M1.obj", Mesh1_V, Mesh1_F);
        saveMeshToFile(shiftMeshAndCurve + "/M2.obj", Mesh2_V, Mesh2_F);
        savePointCloudToFile(shiftMeshAndCurve + "/B1.obj", sortedVertices);
        savePointCloudToFile(shiftMeshAndCurve + "/B2.obj", sortedVertices2);

        polyscope::init();
        polyscope::registerSurfaceMesh("M1", Mesh1_V, Mesh1_F);
        polyscope::registerSurfaceMesh("M2", Mesh2_V, Mesh2_F);
        polyscope::registerPointCloud("B1", sortedVertices);
        polyscope::registerPointCloud("B2", sortedVertices2);
        polyscope::show();


    }

    std::cout << "Step 4 ---------- Lift the 2d boundray curve to 3d []" << std::endl;

    const std::string DEFAULT_2dto3d_FOLDER = "2d_Curve_in_3d";
    if (correspondences2dto3d) {


        const std::string meshesFolder = shiftMeshAndCurve;
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


        auto count = 0;
        auto count2 = 0;
        std::cout << "match V1" << std::endl;

        // Compute the average Z-coordinate for Mesh1
        double averageZ1 = Mesh1_V.col(2).mean(); // Assuming Mesh1_V is Eigen::MatrixXd (Nx3)

        // Set the Z-coordinate of all vertices in V1_pointclouds to averageZ1
        for (Eigen::MatrixXd& V1 : V1_pointclouds) {
            count += V1.rows();
            V1.col(2).setConstant(averageZ1); // Set all Z values of V1 to averageZ1
        }

        // Project and replace vertices for each point cloud in V1_pointclouds
        Eigen::MatrixXd Mesh1_V_new = Mesh1_V; // Start with the original mesh
        Eigen::MatrixXi Mesh1_F_new = Mesh1_F; // Original faces
        std::vector<Eigen::MatrixXd> V1_pointclouds_new; // Store the updated point clouds

        for (Eigen::MatrixXd& V1 : V1_pointclouds) {
            Eigen::MatrixXd V1_pointcloud_new; // Temporary variable for updated point cloud
            projectAndReplaceVertices(Mesh1_V_new, Mesh1_F_new, V1, V1_pointcloud_new, Mesh1_V_new, Mesh1_F_new);
            V1_pointclouds_new.push_back(V1_pointcloud_new); // Store the updated point cloud
            count2 += V1.rows();
        }

        std::cout << "After matching V1, size difference: " << std::to_string(count) << " -- to -- " << std::to_string(count2 )<< std::endl;


        std::cout << "match V2" << std::endl;

        Eigen::MatrixXd Mesh2_V_new = Mesh2_V; // Start with the original mesh
        Eigen::MatrixXi Mesh2_F_new = Mesh2_F;
        std::vector<Eigen::MatrixXd> V2_pointclouds_new;


        double averageZ = Mesh2_V.col(2).mean();
        for (Eigen::MatrixXd& V2 : V2_pointclouds) {
            V2.col(2).setConstant(averageZ); // Set all Z values of V2 to averageZ
        }
        polyscope::init();
        polyscope::registerSurfaceMesh("Before", Mesh2_V, Mesh2_F);
        int i = 0;

        for (Eigen::MatrixXd& V2 : V2_pointclouds) {
            polyscope::registerPointCloud("before_" + std::to_string(i), V2);
            i += 1;
        }


        polyscope::show();


        auto count_fake = 0;
        auto count_fake2 = 0;
        for (Eigen::MatrixXd& V2 : V2_pointclouds) {
            count_fake += V2.rows();

            Eigen::MatrixXd V2_pointcloud_new; // Create a temporary for the new point cloud
            // Project and split mesh, update V2_pointcloud_new, Mesh2_V_new, and Mesh2_F_new
            //projectAndSplitMesh(Mesh2_V_new, Mesh2_F_new, V2, V2_pointcloud_new, Mesh2_V_new, Mesh2_F_new);
            projectAndReplaceVertices(Mesh2_V_new, Mesh2_F_new, V2, V2_pointcloud_new, Mesh2_V_new, Mesh2_F_new);
            // Store the updated point cloud in the vector
            V2_pointclouds_new.push_back(V2_pointcloud_new);

            //findClosestCorrespondences(Mesh2_V, V2);
            //V2_pointclouds_new.push_back(V2);
            count_fake2 += V2.rows();


        }

        std::cout << "After matching V2, size difference: " << std::to_string(count_fake) << " -- to -- " << std::to_string(count_fake2) << std::endl;

        std::cout << "Mesh Old size: " << Mesh2_V.rows() << std::endl;
        std::cout << "Faces Old size: " << Mesh2_F.rows() << std::endl;
        std::cout << "Mesh new size: " << Mesh2_V_new.rows() << std::endl;
        std::cout << "Faces new size: " << Mesh2_F_new.rows() << std::endl;




        // Check if meshes were successfully read
        if (Mesh1_V.rows() == 0 || Mesh2_V.rows() == 0) {
            std::cerr << "Error: Meshes not successfully loaded." << std::endl;
            return -1;
        }

        writeOutputsToFolder(Curve2dTo3dFolder, Mesh1_V_new, Mesh1_F_new, Mesh2_V_new, Mesh2_F_new, V1_pointclouds_new, V2_pointclouds_new);
        //saveMeshToFile(DEFAULT_2dto3d_FOLDER + "/V2_new.obj", Mesh2_V_new, Mesh2_F_new);
 
        // Show them in Polyscope with the common color
        showInPolyscope(Mesh1_V_new, Mesh1_F_new, Mesh2_V_new, Mesh2_F_new, V1_pointclouds_new, V2_pointclouds_new);


        int numHoles = countHoles(Mesh2_V_new, Mesh2_F_new);
        std::cout << "Number of holes in the mesh: " << numHoles << std::endl;


    }

    std::cout << "Step 5 ---------- Parameterize surface using lifted curve as boundary" << std::endl;

    const std::string surfaceParam = "surfaceParameterize";
    if (parameterizeSurfaceBool) {

        std::cout << "surface Parameterization:" << std::endl;

        createDirectory(surfaceParam);
        clearDirectory(surfaceParam);

        Eigen::MatrixXd Mesh1_V;
        Eigen::MatrixXi Mesh1_F;
        readMeshFromFile(DEFAULT_2dto3d_FOLDER + "/Mesh1.obj", Mesh1_V, Mesh1_F);

        Eigen::MatrixXd Mesh2_V;
        Eigen::MatrixXi Mesh2_F;
        readMeshFromFile(DEFAULT_2dto3d_FOLDER + "/Mesh2.obj", Mesh2_V, Mesh2_F);

        std::cout << "Mesh 2 size: " << Mesh2_V.rows() << std::endl;
        std::cout << "Faces 2 size: " << Mesh2_F.rows() << std::endl;


        std::string V1_regex = "V1_pointcloud_(\\d+)\\.txt";
        std::string V2_regex = "V2_pointcloud_(\\d+)\\.txt";

        Eigen::MatrixXd V3 = readAndConcatenatePointClouds(DEFAULT_2dto3d_FOLDER, V1_regex);
        auto xx = V3.rows();
        std::vector<int> duplicatesV1 = std::vector<int>();
        findExactCorrespondences(Mesh1_V, V3, duplicatesV1);
        auto xx2 = V3.rows();


        std::cout << "----------------After matching V1, size difference: " << std::to_string(xx) << " -- to -- " << std::to_string(xx2) <<  " length: " << duplicatesV1.size() << std::endl;


        Eigen::MatrixXd V3_obj2 = readAndConcatenatePointClouds(DEFAULT_2dto3d_FOLDER, V2_regex);

        for (int j = duplicatesV1.size() - 1; j >= 0; j--) {
            int row_to_remove = duplicatesV1[j];
            if (row_to_remove >= 0 && row_to_remove < V3_obj2.rows()) {
                Eigen::MatrixXd temp(V3_obj2.rows() - 1, V3_obj2.cols());
                temp << V3_obj2.topRows(row_to_remove),
                    V3_obj2.bottomRows(V3_obj2.rows() - row_to_remove - 1);
                V3_obj2 = temp;
            }
        }

        std::vector<int> duplicatesV2 = std::vector<int>();
        auto yy = V3_obj2.rows();
        findExactCorrespondences(Mesh2_V, V3_obj2, duplicatesV2);
        auto yy2 = V3_obj2.rows();

        for (int j = duplicatesV2.size() - 1; j >= 0; j--) {
            int row_to_remove = duplicatesV2[j];
            if (row_to_remove >= 0 && row_to_remove < V3.rows()) {
                Eigen::MatrixXd temp(V3.rows() - 1, V3.cols());
                temp << V3.topRows(row_to_remove),
                    V3.bottomRows(V3.rows() - row_to_remove - 1);
                V3 = temp;
            }
        }


        if (true) {
            std::cout << "surface parameterization (MESH 1) using LCSM without splititng the mesh into 2" << std::endl;
            Eigen::MatrixXd UV_map;
            polyscope::init();
            polyscope::options::programName = "No Split Mesh LCSM, projection";
            polyscope::registerPointCloud("borderV3", V3);
            if (!paramsurface5(Mesh1_V, Mesh1_F, UV_map, V3, true, V3)) {
                std::cerr << "Surface parameterization failed.\n";
                return EXIT_FAILURE;
            }

            saveMeshToFile(surfaceParam + "/M1.obj", Mesh1_V, Mesh1_F);
            write2DVerticesToPLY(surfaceParam + "/UV1.obj", UV_map);
            writeVerticesToPLY(surfaceParam + "/B1.obj", V3);
        }

        if (true) {
            std::cout << "surface parameterization (MESH 2) using LCSM without splititng the mesh into 2" << std::endl;
            Eigen::MatrixXd UV_map;
            polyscope::init();

            polyscope::init();
            polyscope::options::programName = "No Split Mesh LCSM, projection";
            polyscope::registerPointCloud("border V3 shift", V3_obj2);


            CalculateGenus(Mesh2_V, Mesh2_F);




            if (!paramsurface5(Mesh2_V, Mesh2_F, UV_map, V3_obj2, true, V3)) {
                std::cerr << "Surface parameterization failed.\n";
                return EXIT_FAILURE;
            }

            saveMeshToFile(surfaceParam + "/M2.obj", Mesh2_V, Mesh2_F);
            write2DVerticesToPLY(surfaceParam + "/UV2.obj", UV_map);
            writeVerticesToPLY(surfaceParam + "/B2.obj", V3_obj2);


        }

  

        if (false) {

            std::cout << "split Meshes" << std::endl;

            const std::string splitmesh = "SplitMeshes";
            createDirectory(splitmesh);
            clearDirectory(splitmesh);


            Eigen::MatrixXd MeshA_V;
            Eigen::MatrixXi MeshA_F;
            Eigen::MatrixXd Mesh1_A_Removed_V;
            Eigen::MatrixXi Mesh1_A_Removed_F;

            
            std::cout << "Step 1: remove lifted curve vertices and faces" << std::endl;
            RemoveVerticesAndFaces(Mesh1_V, Mesh1_F, V3, MeshA_V, MeshA_F, Mesh1_A_Removed_V, Mesh1_A_Removed_F);
            countConnectedComponents(MeshA_F);

            Eigen::MatrixXd border_V = getBorderVerticesMatrix(MeshA_V, MeshA_F);

            //showMeshAndPointCloud(MeshA_V, MeshA_F, border_V);


            Eigen::MatrixXd MeshB_V;
            Eigen::MatrixXi MeshB_F;
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
            //FindMatchingEdges(MeshA_V, MeshA_F, removed_V, removed_F, 1e-4, Mesh_V_Split1, Mesh_F_Split1);
            //FindMatchingEdges(MeshA_V, MeshA_F, removed_V, removed_F, 1e-4, Mesh_V_Split1, Mesh_F_Split1);
            //FindMatchingEdges(MeshA_V, MeshA_F, removed_V, removed_F, 1e-4, Mesh_V_Split2, Mesh_F_Split2);
            //FindMatchingEdges(MeshA_V, MeshA_F, removed_V, removed_F, 1e-4, Mesh_V_Split2, Mesh_F_Split2);

            //polyscope::init();
            //polyscope::registerSurfaceMesh("Mesh INput", Mesh1_V, Mesh1_F);
            //polyscope::registerSurfaceMesh("Mesh AAA", MeshA_V, MeshA_F);
            //polyscope::registerSurfaceMesh("Mesh BBB", MeshB_V, MeshB_F);
            //polyscope::registerSurfaceMesh("Mesh P1", Mesh_V_Split1, Mesh_F_Split1);
            //polyscope::registerSurfaceMesh("Mesh P2", Mesh_V_Split2, Mesh_F_Split2);
            //polyscope::registerPointCloud("INput lifted curve", V3);
            //polyscope::registerPointCloud("Removed A to B", removed_V);
            //polyscope::registerPointCloud("Border A to B", border_V_new);
            //polyscope::show();


            //saveMeshToFile(splitmesh + "/A.obj", Mesh_V_Split1, Mesh_F_Split1);
            //saveMeshToFile(splitmesh + "/B.obj", Mesh_V_Split2, Mesh_F_Split2);


            //// parameterize whole surface
            if (true) {
                Eigen::MatrixXd UV_split1;
                Eigen::MatrixXd V_border_split1 = getBorderVerticesMatrix(Mesh_V_Split1, Mesh_F_Split1);
                polyscope::options::programName = "Split Mesh (one side) LCSM, projection";
                if (!paramsurface5(Mesh_V_Split1, Mesh_F_Split1, UV_split1, V_border_split1, true, V_border_split1)) {
                    std::cerr << "Surface parameterization failed.\n";
                    return EXIT_FAILURE;
                }
            }


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


    }


    const std::string correspondence3dMatched = "correspondence3dMatched";

    if (uvMapCorrespondence) {


        Eigen::MatrixXd V1; // Mesh 1 vertices       n x 3
        Eigen::MatrixXi F1; // Mesh 1 faces          m x 3
        Eigen::MatrixXd B1; // boundary 1 vertices   p x 3
        Eigen::MatrixXd UV1; // uv map of mesh 1     n x 2

        Eigen::MatrixXd V2; // mesh 2 vertices       s x 3
        Eigen::MatrixXi F2; // mesh 2 faces          q x 3
        Eigen::MatrixXd B2; // boundary 2 vertices   p x 3
        Eigen::MatrixXd UV2; // uv map of mesh 2     s x 2

        readMeshFromFile(surfaceParam + "/M1.obj", V1, F1);
        readMeshFromFile(surfaceParam + "/M2.obj", V2, F2);
        B1 = readVerticesFromPLY(surfaceParam + "/B1.obj");
        B2 = readVerticesFromPLY(surfaceParam + "/B2.obj");
        UV1 = readVerticesFromPLY2D(surfaceParam + "/UV1.obj");
        UV2 = readVerticesFromPLY2D(surfaceParam + "/UV2.obj");

        // Debug prints
        std::cout << "V1: " << V1.rows() << " x " << V1.cols() << std::endl;
        std::cout << "F1: " << F1.rows() << " x " << F1.cols() << std::endl;
        std::cout << "B1: " << B1.rows() << " x " << B1.cols() << std::endl;
        std::cout << "UV1: " << UV1.rows() << " x " << UV1.cols() << std::endl;

        std::cout << "V2: " << V2.rows() << " x " << V2.cols() << std::endl;
        std::cout << "F2: " << F2.rows() << " x " << F2.cols() << std::endl;
        std::cout << "B2: " << B2.rows() << " x " << B2.cols() << std::endl;
        std::cout << "UV2: " << UV2.rows() << " x " << UV2.cols() << std::endl;


        createDirectory(correspondence3dMatched);
        UVToCorrespondence(V1, F1, B1, UV1, V2, F2, B2, UV2, correspondence3dMatched);
    }

    const std::string evaluateCorrespondenceFolder = "EvaluteFolder";

    if (evaluateCorrespondence) {

        

        Eigen::MatrixXd V1; // Mesh 1 vertices       n x 3
        Eigen::MatrixXi F1; // Mesh 1 faces          m x 3
        Eigen::MatrixXd V2; // mesh 2 vertices       s x 3
        Eigen::MatrixXi F2; // mesh 2 faces          q x 3
        readMeshFromFile(correspondence3dMatched + "/M1.obj", V1, F1);
        readMeshFromFile(correspondence3dMatched + "/M2.obj", V2, F2);


        analyzeAndVisualizeCorrespondence(V1, F1, V2, F2, evaluateCorrespondenceFolder);

    }



    if (newShapeMake) {
        Eigen::MatrixXd V1, V2, V3;
        Eigen::MatrixXi F1, F2, F3;

        readMeshFromFile(correspondence3dMatched + "/M1.obj", V1, F1);
        readMeshFromFile(correspondence3dMatched + "/M2.obj", V2, F2);
        readMeshFromFile(correspondence3dMatched + "/M2.obj", V3, F3);

        // Combine into a vector of pairs
        std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> inputShapes = {
            {V1, F1},
            {V2, F2},
            {V3, F3}
        };

        // Call the main_phase2 function
        main_phase2(inputShapes);

    }



    return 0;
}



