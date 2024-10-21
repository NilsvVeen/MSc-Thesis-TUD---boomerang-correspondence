// parameterizeSurface.cpp

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/boost/graph/properties.h>
#include <iostream>
#include <fstream>
#include "libraries/CGAL-5.6.1/include/CGAL/Polygon_mesh_processing/measure.h"

// Typedefs
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor;

int main(int argc, char* argv[]) {
    // Ensure file path is provided
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " input_mesh.obj\n";
        return EXIT_FAILURE;
    }

    // Create surface mesh and read from file
    Surface_mesh mesh;
    std::ifstream input(argv[1]);
    if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
        std::cerr << "Invalid mesh file.\n";
        return EXIT_FAILURE;
    }

    // Ensure the mesh has a border (necessary for parameterization)
    if (CGAL::Polygon_mesh_processing::longest_border(mesh).first == Surface_mesh::null_halfedge()) {
        std::cerr << "Input mesh has no border.\n";
        return EXIT_FAILURE;
    }

    // Property map for UV coordinates
    Surface_mesh::Property_map<Surface_mesh::Vertex_index, Kernel::Point_2> uvmap =
        mesh.add_property_map<Surface_mesh::Vertex_index, Kernel::Point_2>("v:uv").first;

    // Get a border halfedge
    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh).first;

    // Define parameterizer types
    typedef CGAL::Surface_mesh_parameterization::Circular_border_parameterizer_3<Surface_mesh> Border_parameterizer;
    typedef CGAL::Surface_mesh_parameterization::Discrete_authalic_parameterizer_3<Surface_mesh, Border_parameterizer> Parameterizer;

    // Perform the parameterization
    CGAL::Surface_mesh_parameterization::Error_code status =
        CGAL::Surface_mesh_parameterization::parameterize(mesh, bhd, uvmap);

    // Check for success
    if (status != CGAL::Surface_mesh_parameterization::OK) {
        std::cerr << "Parameterization failed: " << status << std::endl;
        return EXIT_FAILURE;
    }

    // Output UV coordinates to file
    std::ofstream uv_output("parameterized_mesh.obj");
    for (auto v : mesh.vertices()) {
        const Kernel::Point_2& uv = uvmap[v];
        uv_output << "vt " << uv.x() << " " << uv.y() << std::endl;
    }

    std::cout << "Parameterization successful, UV coordinates saved to 'parameterized_mesh.obj'.\n";
    return EXIT_SUCCESS;
}
