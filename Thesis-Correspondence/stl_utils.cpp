#include <iostream>
#include <filesystem>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Dense>

//#include <opencv2/opencv.hpp> // For saving the PNG

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>

#include <igl/readSTL.h>
#include <igl/writeSTL.h>
#include <igl/fit_plane.h>
#include <igl/slice_mask.h> // To filter the vertices
#include <igl/copyleft/cgal/convex_hull.h>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef K::FT                                                FT;
typedef K::Point_2                                           Point;
typedef K::Segment_2                                         Segment;
typedef CGAL::Alpha_shape_vertex_base_2<K>                   Vb;
typedef CGAL::Alpha_shape_face_base_2<K>                     Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>          Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>                 Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator            Alpha_shape_edges_iterator;





void viewSTLObject(const std::string& filename) {
    Eigen::MatrixXd V; // Vertices
    Eigen::MatrixXi F; // Faces
    Eigen::MatrixXd N; // Normals (if required)

    // Open the STL file using ifstream
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Failed to open STL file: " << filename << std::endl;
        return;
    }

    // Read the STL file using libigl
    if (!igl::readSTL(file, V, F, N)) {
        std::cerr << "Failed to load STL file: " << filename << std::endl;
        return;
    }

    // Close the file
    file.close();

    // Initialize the Polyscope viewer
    polyscope::init();

    // Register the mesh to Polyscope
    polyscope::registerSurfaceMesh("STL Mesh", V, F);

    // Show the viewer
    polyscope::show();
}

void printFilesInDirectory(const std::string& directoryPath) {
    try {
        std::filesystem::path dir(directoryPath);
        if (!std::filesystem::exists(dir)) {
            std::cerr << directoryPath << " does not exist." << std::endl;
            return;
        }
        if (!std::filesystem::is_directory(dir)) {
            std::cerr << directoryPath << " is not a valid directory." << std::endl;
            return;
        }

        std::cout << "Files in directory " << directoryPath << ":" << std::endl;
        for (const auto& entry : std::filesystem::directory_iterator(dir)) {
            std::cout << entry.path().filename().string() << std::endl;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error accessing directory: " << e.what() << std::endl;
    }
}


// Save mesh to a file in OBJ format
void saveMeshToFile(const std::string& filename, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    // Write vertices
    for (int i = 0; i < V.rows(); ++i) {
        file << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
    }

    // Write faces
    for (int i = 0; i < F.rows(); ++i) {
        file << "f " << F(i, 0) + 1 << " " << F(i, 1) + 1 << " " << F(i, 2) + 1 << std::endl;
    }

    file.close();
    std::cout << "Mesh saved to " << filename << std::endl;
}

// Save point cloud to a file in PLY format
void savePointCloudToFile(const std::string& filename, const Eigen::MatrixXd& V) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    // Write PLY header
    file << "ply\nformat ascii 1.0\nelement vertex " << V.rows() << "\n";
    file << "property float x\nproperty float y\nproperty float z\nend_header\n";

    // Write vertices
    for (int i = 0; i < V.rows(); ++i) {
        file << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
    }

    file.close();
    std::cout << "Point cloud saved to " << filename << std::endl;
}



// Function to compute the surface area of a triangle
double computeTriangleArea(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) {
    return 0.5 * (v1 - v0).cross(v2 - v0).norm();
}

// Project vertices to 2D (XY plane) and fix Z to a specified value
Eigen::MatrixXd projectToXYWithFixedZ(const Eigen::MatrixXd& V, double fixedZ) {
    Eigen::MatrixXd V_2D(V.rows(), 3);
    V_2D.leftCols(2) = V.leftCols(2);  // Copy X and Y
    V_2D.col(2).setConstant(fixedZ);  // Set Z to fixedZ
    return V_2D;
}

// Initialize Polyscope and register a mesh
void initializePolyscopeAndRegisterMesh(const std::string& name, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    polyscope::init();
    polyscope::registerSurfaceMesh(name, V, F);
}

// Take a screenshot of the current Polyscope view
void takeScreenshot(const std::string& screenshotName) {
    polyscope::screenshot(screenshotName);
}

// Show all registered meshes in Polyscope
void showMeshes() {
    polyscope::show();
}

// Function to fit a plane to the mesh and align it
void fitPlaneAndAlignMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& rotatedV, Eigen::MatrixXi& rotatedF) {
    Eigen::RowVector3d planeNormal;
    Eigen::RowVector3d planePoint;
    igl::fit_plane(V, planeNormal, planePoint);

    Eigen::Vector3d normal(planeNormal);
    Eigen::Vector3d targetDirection(0, 0, 1);
    Eigen::Matrix3d rotationMatrix = Eigen::Quaterniond::FromTwoVectors(normal, targetDirection).toRotationMatrix();

    rotatedV = (rotationMatrix * V.transpose()).transpose();
    rotatedF = F;
}

// Function to find the slice with maximum area
void findMaxAreaSlice(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& V_maxAreaSlice, Eigen::MatrixXi& F_maxAreaSlice, int& maxAreaSlice) {
    int numSlices = 10;
    double minZ = V.col(2).minCoeff();
    double maxZ = V.col(2).maxCoeff();
    double sliceThickness = (maxZ - minZ) / numSlices;

    double maxSliceArea = 0.0;
    maxAreaSlice = -1;

    for (int i = 0; i < numSlices; ++i) {
        double sliceMinZ = minZ + i * sliceThickness;
        double sliceMaxZ = sliceMinZ + sliceThickness;

        std::vector<int> sliceFaces;

        for (int j = 0; j < F.rows(); ++j) {
            Eigen::Vector3d v0 = V.row(F(j, 0));
            Eigen::Vector3d v1 = V.row(F(j, 1));
            Eigen::Vector3d v2 = V.row(F(j, 2));

            if ((v0.z() >= sliceMinZ && v0.z() < sliceMaxZ) &&
                (v1.z() >= sliceMinZ && v1.z() < sliceMaxZ) &&
                (v2.z() >= sliceMinZ && v2.z() < sliceMaxZ)) {
                sliceFaces.push_back(j);
            }
        }

        Eigen::MatrixXi F_slice(sliceFaces.size(), 3);
        for (int k = 0; k < sliceFaces.size(); ++k) {
            F_slice.row(k) = F.row(sliceFaces[k]);
        }

        double sliceArea = 0.0;
        for (int j = 0; j < F_slice.rows(); ++j) {
            Eigen::Vector3d v0 = V.row(F_slice(j, 0));
            Eigen::Vector3d v1 = V.row(F_slice(j, 1));
            Eigen::Vector3d v2 = V.row(F_slice(j, 2));
            sliceArea += computeTriangleArea(v0, v1, v2);
        }

        if (sliceArea > maxSliceArea) {
            maxSliceArea = sliceArea;
            maxAreaSlice = i;
            V_maxAreaSlice = V;
            F_maxAreaSlice = F_slice;
        }
    }
}

// Function to calculate average Z coordinate for a slice
double calculateAverageZ(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    Eigen::VectorXd zCoords(F.rows());
    for (int i = 0; i < F.rows(); ++i) {
        zCoords(i) = V.row(F(i, 0)).z() +
            V.row(F(i, 1)).z() +
            V.row(F(i, 2)).z();
    }
    return zCoords.mean() / 3.0;
}

// Function to handle STL file loading and processing
void processSTLFile(const std::string& filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& N) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Failed to open STL file: " << filename << std::endl;
        return;
    }

    if (!igl::readSTL(file, V, F, N)) {
        std::cerr << "Failed to load STL file: " << filename << std::endl;
        return;
    }

    file.close();
}


void findBorderVerticesWithAlphaShape(const Eigen::MatrixXd& V_2D, std::vector<Eigen::Vector2d>& borderVertices, double& alpha) {
    std::vector<Point> points;
    for (int i = 0; i < V_2D.rows(); ++i) {
        points.emplace_back(V_2D(i, 0), V_2D(i, 1));
    }

    // Create the Alpha Shape with an initial alpha parameter
    Alpha_shape_2 alpha_shape(points.begin(), points.end(), FT(0.0), Alpha_shape_2::GENERAL);
    
    // Set an appropriate alpha value for better border extraction (adjust as needed)
    alpha_shape.set_alpha(alpha); // Modify the alpha parameter based on your data for the best results

    // Extract border edges
    std::set<Point> borderPoints;
    for (Alpha_shape_edges_iterator it = alpha_shape.alpha_shape_edges_begin(); it != alpha_shape.alpha_shape_edges_end(); ++it) {
        if (alpha_shape.classify(*it) == Alpha_shape_2::REGULAR || alpha_shape.classify(*it) == Alpha_shape_2::SINGULAR) {
            Segment segment = alpha_shape.segment(*it);
            borderPoints.insert(segment.source());
            borderPoints.insert(segment.target());
        }
    }

    // Convert border points back to Eigen format
    for (const auto& point : borderPoints) {
        borderVertices.emplace_back(point.x(), point.y());
    }
}

// Main function to fit plane, align mesh, and show results
std::vector<Eigen::Vector2d> fitPlaneAndAlignMesh(const std::string& filename, const std::string& outputDir) {
    Eigen::MatrixXd V; // Vertices
    Eigen::MatrixXi F; // Faces
    Eigen::MatrixXd N; // Normals
    std::vector<Eigen::Vector2d> borderVertices;

    // Process STL file
    processSTLFile(filename, V, F, N);

    // Fit plane and align mesh
    Eigen::MatrixXd rotatedV;
    Eigen::MatrixXi rotatedF;
    fitPlaneAndAlignMesh(V, F, rotatedV, rotatedF);



    // Define file paths for saving the meshes
    std::string rotatedMeshFile = outputDir + "/rotated_mesh.obj";
    std::string maxAreaSliceFile = outputDir + "/max_area_slice.obj";
    std::string projectionFile = outputDir + "/2d_projection.obj";
    std::string borderVerticesFile = outputDir + "/alpha_shape_border.obj";

    // Save the rotated mesh to a file
    saveMeshToFile(rotatedMeshFile, rotatedV, rotatedF);

    // Initialize Polyscope and register the rotated mesh
    initializePolyscopeAndRegisterMesh("Rotated Mesh", rotatedV, rotatedF);
    takeScreenshot("rotated_mesh.png");

    // Find the slice with maximum area
    Eigen::MatrixXd V_maxAreaSlice;
    Eigen::MatrixXi F_maxAreaSlice;
    int maxAreaSlice;
    findMaxAreaSlice(rotatedV, rotatedF, V_maxAreaSlice, F_maxAreaSlice, maxAreaSlice);

    if (maxAreaSlice >= 0) {
        // Calculate average Z coordinate for the max area slice
        double avgZ = calculateAverageZ(rotatedV, F_maxAreaSlice);

        // Project the mesh onto the XY plane with fixed Z coordinate
        Eigen::MatrixXd V_2D = projectToXYWithFixedZ(rotatedV, avgZ);

        // Register the 2D projection and max area slice with Polyscope
        initializePolyscopeAndRegisterMesh("2D Projection", V_2D, rotatedF);
        initializePolyscopeAndRegisterMesh("Max Area Slice", V_maxAreaSlice, F_maxAreaSlice);

        // Save the 2D projection and max area slice to files
        saveMeshToFile(maxAreaSliceFile, V_maxAreaSlice, F_maxAreaSlice);
        saveMeshToFile(projectionFile, V_2D, rotatedF);

        // Create a vector for the convex hull vertices
        Eigen::MatrixXd V_hull;
        Eigen::MatrixXi F_hull; // Faces of the convex hull

        // Call the convex hull function from libigl
        igl::copyleft::cgal::convex_hull(V_2D, V_hull, F_hull);

        // Register the convex hull with Polyscope
        initializePolyscopeAndRegisterMesh("Convex Hull", V_hull, F_hull);

        // Try 2D Alpha Shape from CGAL
        double alpha = 1.0; // Alpha for the alpha shape
        findBorderVerticesWithAlphaShape(V_2D, borderVertices, alpha);

        if (!borderVertices.empty()) {
            // Convert border vertices to Eigen format and set last column to avgZ
            Eigen::MatrixXd V_border(borderVertices.size(), 3);
            for (size_t i = 0; i < borderVertices.size(); ++i) {
                V_border(i, 0) = borderVertices[i].x();
                V_border(i, 1) = borderVertices[i].y();
                V_border(i, 2) = avgZ;
            }

            // Register the border vertices as a point cloud with Polyscope
            std::string cloudName = "Alpha Shape Border Vertices (alpha = " + std::to_string(alpha) + ")";
            polyscope::registerPointCloud(cloudName, V_border);
            takeScreenshot("alpha_shape_border_alpha_" + std::to_string(alpha) + ".png");

            // Save the border vertices to a file
            savePointCloudToFile(borderVerticesFile, V_border);
        }

        // Take screenshots
        takeScreenshot("2d_projection.png");
        takeScreenshot("max_area_slice.png");
    }

    // Show the registered meshes
    showMeshes();

    return borderVertices;
}
