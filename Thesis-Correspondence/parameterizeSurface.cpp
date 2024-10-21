#include "parameterizeSurface.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/boost/graph/properties.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Surface_mesh_parameterization/LSCM_parameterizer_3.h>

// Typedefs
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor;

void convertEigenToSurfaceMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Surface_mesh& mesh) {
    std::vector<Surface_mesh::Vertex_index> vertices;

    // Add vertices
    for (int i = 0; i < V.rows(); ++i) {
        vertices.push_back(mesh.add_vertex(Kernel::Point_3(V(i, 0), V(i, 1), V(i, 2))));
        //std::cout << "Added vertex " << i << ": (" << V(i, 0) << ", " << V(i, 1) << ", " << V(i, 2) << ")\n";
    }

    // Add faces
    for (int i = 0; i < F.rows(); ++i) {
        int v1 = F(i, 0), v2 = F(i, 1), v3 = F(i, 2);
        if (v1 >= vertices.size() || v2 >= vertices.size() || v3 >= vertices.size()) {
            std::cerr << "Invalid face index at row " << i << ": (" << v1 << ", " << v2 << ", " << v3 << ")\n";
            return; // or handle the error
        }
        mesh.add_face(vertices[v1], vertices[v2], vertices[v3]);
        //std::cout << "Added face " << i << ": (" << v1 << ", " << v2 << ", " << v3 << ")\n";
    }
}



// Function to parameterize surface using LSCM
bool parameterizeSurface(const Eigen::MatrixXd& Mesh1_V, const Eigen::MatrixXi& Mesh1_F, const std::string& output_filename) {
    // Create CGAL surface mesh
    Surface_mesh mesh;
    convertEigenToSurfaceMesh(Mesh1_V, Mesh1_F, mesh);

    // Ensure the mesh has a border
    if (CGAL::Polygon_mesh_processing::longest_border(mesh).first == Surface_mesh::null_halfedge()) {
        std::cerr << "Mesh has no border.\n";
        return false;
    }

    // Property map for UV coordinates
    Surface_mesh::Property_map<Surface_mesh::Vertex_index, Kernel::Point_2> uvmap =
        mesh.add_property_map<Surface_mesh::Vertex_index, Kernel::Point_2>("v:uv").first;

    // Property map for vertex indices
    Surface_mesh::Property_map<Surface_mesh::Vertex_index, int> vertex_index_map =
        mesh.add_property_map<Surface_mesh::Vertex_index, int>("v:index").first;

    // Property map for parameterized vertices (not strictly necessary but can be useful)
    Surface_mesh::Property_map<Surface_mesh::Vertex_index, bool> vertex_parameterized_map =
        mesh.add_property_map<Surface_mesh::Vertex_index, bool>("v:parameterized").first;

    // Get a border halfedge
    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh).first;

    // Define parameterizer type for LSCM
    using LSCM_Parameterizer = CGAL::Surface_mesh_parameterization::LSCM_parameterizer_3<Surface_mesh>;

    // Create the LSCM parameterizer instance
    LSCM_Parameterizer parameterizer;

    // Perform the LSCM parameterization
    CGAL::Surface_mesh_parameterization::Error_code status =
        parameterizer.parameterize(mesh, bhd, uvmap, vertex_index_map, vertex_parameterized_map);

    // Check for success
    if (status != CGAL::Surface_mesh_parameterization::OK) {
        std::cerr << "Parameterization failed: " << status << std::endl;
        return false;
    }

    // Output UV coordinates to file
    std::ofstream uv_output(output_filename);
    for (auto v : mesh.vertices()) {
        const Kernel::Point_2& uv = uvmap[v];
        uv_output << "vt " << uv.x() << " " << uv.y() << std::endl;
    }

    std::cout << "Parameterization successful, UV coordinates saved to '" << output_filename << "'.\n";
    return true;
}





