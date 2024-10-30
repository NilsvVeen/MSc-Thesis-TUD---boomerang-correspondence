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