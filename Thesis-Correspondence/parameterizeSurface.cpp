#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_parameterization/LSCM_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>          Kernel;
typedef Kernel::Point_2                         Point_2;
typedef Kernel::Point_3                         Point_3;
typedef CGAL::Surface_mesh<Point_3>             SurfaceMesh;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor     vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor   halfedge_descriptor;

namespace SMP = CGAL::Surface_mesh_parameterization;

//bool parameterizeSurface(const Eigen::MatrixXd& Mesh1_V, const Eigen::MatrixXi& Mesh1_F, const std::string& output_filename) {
//    SurfaceMesh sm;
//
//    // Convert Eigen matrices to CGAL surface mesh
//    for (int i = 0; i < Mesh1_V.rows(); ++i) {
//        sm.add_vertex(Point_3(Mesh1_V(i, 0), Mesh1_V(i, 1), Mesh1_V(i, 2)));
//    }
//    for (int i = 0; i < Mesh1_F.rows(); ++i) {
//        std::vector<vertex_descriptor> vertices;
//        for (int j = 0; j < Mesh1_F.cols(); ++j) {
//            vertices.push_back(vertex_descriptor(Mesh1_F(i, j)));
//        }
//        sm.add_face(vertices);
//    }
//
//    // Verify the mesh is valid
//    if (!CGAL::is_valid(sm)) {
//        std::cerr << "Input mesh is invalid." << std::endl;
//        return false;
//    }
//
//    // Get a halfedge on the border
//    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;
//
//    // Check if bhd is valid
//    if (!sm.is_border(bhd)) {
//        std::cerr << "The specified halfedge is not a border halfedge." << std::endl;
//        return false;
//    }
//
//    // Iterate through the boundary loop to find two vertices
//    std::vector<vertex_descriptor> boundary_vertices;
//    halfedge_descriptor hd = bhd;
//    do {
//        boundary_vertices.push_back(sm.target(hd));
//        hd = sm.next(hd);
//    } while (hd != bhd && boundary_vertices.size() < 2); // Stop after finding two vertices
//
//    if (boundary_vertices.size() < 2) {
//        std::cerr << "Failed to find two distinct boundary vertices." << std::endl;
//        return false;
//    }
//
//    // Create UV property map for the parameterized values
//    typedef SurfaceMesh::Property_map<vertex_descriptor, Point_2> UV_pmap;
//    UV_pmap uv_map = sm.add_property_map<vertex_descriptor, Point_2>("h:uv").first;
//
//    // Create a vertex index property map
//    typedef SurfaceMesh::Property_map<vertex_descriptor, std::size_t> Index_map;
//    Index_map index_map = sm.add_property_map<vertex_descriptor, std::size_t>("h:index").first;
//    std::size_t index = 0;
//    for (auto vd : sm.vertices()) {
//        index_map[vd] = index++;
//    }
//
//    // Create a vertex parameterized property map
//    typedef SurfaceMesh::Property_map<vertex_descriptor, bool> Parameterized_map;
//    Parameterized_map param_map = sm.add_property_map<vertex_descriptor, bool>("h:parameterized").first;
//
//    // Pin two boundary vertices with fixed UV coordinates
//    auto v1 = boundary_vertices[0];  // First pinned vertex
//    auto v2 = boundary_vertices[1];  // Second pinned vertex
//
//    uv_map[v1] = Point_2(0, 0);  // Pin first vertex to (0, 0)
//    uv_map[v2] = Point_2(1, 0);  // Pin second vertex to (1, 0)
//
//    param_map[v1] = true; // Mark v1 as parameterized (pinned)
//    param_map[v2] = true; // Mark v2 as parameterized (pinned)
//
//    // Create an LSCM parameterizer instance
//    SMP::LSCM_parameterizer_3<SurfaceMesh> lscm_param;
//
//    // Perform the parameterization
//    bool success = lscm_param.parameterize(sm, bhd, uv_map, index_map, param_map);
//    if (!success) {
//        std::cerr << "Parameterization failed." << std::endl;
//        return false; // Parameterization failed
//    }
//
//    // Check UV values for validity
//    for (const auto& vd : sm.vertices()) {
//        const Point_2& uv = uv_map[vd];
//        if (std::isnan(uv.x()) || std::isnan(uv.y())) {
//            std::cerr << "UV for vertex " << vd << " is invalid: (" << uv.x() << ", " << uv.y() << ")\n";
//            return false; // Invalid UV coordinate
//        }
//        std::cout << "UV for vertex " << vd << ": (" << uv.x() << ", " << uv.y() << ")\n";
//    }
//
//    // Write the UV map to the output file
//    std::ofstream out(output_filename);
//    try {
//        SMP::IO::output_uvmap_to_off(sm, bhd, uv_map, out);
//    }
//    catch (const std::exception& e) {
//        std::cerr << "Failed to write output UV map: " << e.what() << std::endl;
//        return false; // Output failed
//    }
//
//    return true; // Success
//}


#include <igl/cotmatrix.h>
#include <igl/vector_area_matrix.h>
#include <igl/repdiag.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>

//bool parameterizeSurfaceLibiGL(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& UV) {
//    // Step 1: Compute cotangent Laplacian matrix
//    Eigen::SparseMatrix<double> L;
//    igl::cotmatrix(V, F, L);
//
//    // Step 2: Compute area matrix
//    Eigen::SparseMatrix<double> A;
//    igl::vector_area_matrix(F, A);
//
//    // Step 3: Replicate the Laplacian matrix for both u and v coordinates
//    Eigen::SparseMatrix<double> L_flat;
//    igl::repdiag(L, 2, L_flat);
//
//    // Step 4: Extend the area matrix similarly for [u, v] coordinates
//    Eigen::SparseMatrix<double> A_flat;
//    A_flat = 2 * A; // Scaling by 2 as per the LSCM formula
//
//    // Final LSCM energy matrix: L_flat - 2 * A
//    Eigen::SparseMatrix<double> Energy = L_flat - A_flat;
//
//    // Step 5: Fix two arbitrary vertices to remove null space
//    Eigen::VectorXi boundary;
//    igl::boundary_loop(F, boundary);
//    if (boundary.size() < 2) {
//        std::cerr << "Mesh must have at least two boundary vertices." << std::endl;
//        return false;
//    }
//
//    // Fix two boundary vertices (this can be any two distinct vertices)
//    int v1 = boundary(0);
//    int v2 = boundary(1);
//
//    // Step 6: Prepare the RHS (set the UV coordinates for pinned vertices)
//    Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(2 * V.rows(), 2);
//    rhs(2 * v1, 0) = 0.0; rhs(2 * v1 + 1, 1) = 0.0; // Pin (u, v) of v1 to (0, 0)
//    rhs(2 * v2, 0) = 1.0; rhs(2 * v2 + 1, 1) = 0.0; // Pin (u, v) of v2 to (1, 0)
//
//    // Step 7: Solve the linear system
//    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
//    solver.compute(Energy);
//
//    if (solver.info() != Eigen::Success) {
//        std::cerr << "Decomposition failed." << std::endl;
//        return false;
//    }
//
//    Eigen::MatrixXd UV_flat = solver.solve(rhs);
//
//    if (solver.info() != Eigen::Success) {
//        std::cerr << "Solving failed." << std::endl;
//        return false;
//    }
//
//    // Step 8: Extract UV coordinates from the flattened solution
//    UV.resize(V.rows(), 2);
//    for (int i = 0; i < V.rows(); ++i) {
//        UV(i, 0) = UV_flat(2 * i, 0);
//        UV(i, 1) = UV_flat(2 * i + 1, 1);
//    }
//
//    return true;
//}


//#include <igl/map_vertices_to_circle.h>
//#include <igl/harmonic.h>
//#include <igl/boundary_loop.h>
//#include <iostream>
//#include <Eigen/Dense>
//
//bool parameterizeSurfaceHarmonic(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& UV) {
//    // Step 1: Get the boundary vertices
//    Eigen::VectorXi boundary;
//    igl::boundary_loop(F, boundary);
//
//    if (boundary.size() == 0) {
//        std::cerr << "Mesh has no boundary, unable to perform harmonic parameterization." << std::endl;
//        return false;
//    }
//
//    // Step 2: Fix boundary vertices to a circle or predefined shape
//    Eigen::MatrixXd boundary_uv(boundary.size(), 2);
//    for (int i = 0; i < boundary.size(); ++i) {
//        double theta = 2.0 * 3.14159265359 * (double)i / (double)boundary.size();
//        boundary_uv(i, 0) = std::cos(theta); // u-coordinate
//        boundary_uv(i, 1) = std::sin(theta); // v-coordinate
//    }
//
//    // Step 3: Perform harmonic parameterization with fixed boundary
//    igl::harmonic(V, F, boundary, boundary_uv, 1, UV);
//
//    return true;
//}


#include <igl/readOFF.h>
//#include <igl/opengl/glfw/Viewer.h>
#include <igl/lscm.h>
#include <iostream>

// Function to compute UV parametrization and visualize the result
bool paramsurface5(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& UV)
{
    // Perform LSCM parametrization
    igl::lscm(V, F, UV);

    // Check if the UV matrix is valid (it shouldn't be empty)
    if (UV.size() == 0)
    {
        std::cerr << "Failed to compute UV parametrization!" << std::endl;
        return false;
    }

    //// Set up the viewer
    //igl::opengl::glfw::Viewer viewer;

    //// Callback for key press events
    //viewer.callback_key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
    //{
    //    if (key == '1')
    //    {
    //        // Plot the 3D mesh
    //        viewer.data().set_mesh(V, F);
    //        viewer.core().align_camera_center(V, F);
    //    }
    //    else if (key == '2')
    //    {
    //        // Plot the mesh in 2D using the UV coordinates as vertex coordinates
    //        viewer.data().set_mesh(UV, F);
    //        viewer.core().align_camera_center(UV, F);
    //    }

    //    viewer.data().compute_normals();
    //    return false;
    //};

    //// Initial mesh view (3D)
    //viewer.data().set_mesh(V, F);
    //viewer.data().set_uv(UV);

    //// Disable wireframe
    //viewer.data().show_lines = false;

    //// Draw checkerboard texture
    //viewer.data().show_texture = true;

    //// Launch the viewer
    //viewer.launch();

    return true;
}

