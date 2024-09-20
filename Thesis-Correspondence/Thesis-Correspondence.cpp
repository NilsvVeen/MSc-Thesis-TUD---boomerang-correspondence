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
            showPointCloudInParts(sortedVertices);
        }
    }

    return 0;
}



void old_1 (std::string& directoryName) {
    //// Define parameters using integers
    int circleRadius = 10;       // Radius of the circle in units (adjusted for text size)
    int numVertices = 100;       // Number of vertices on the circle and ellipse
    int axisRatio = 2;           // Ratio a/b for the ellipse (integer value)

    //// Generate circle vertices
    std::vector<Point> circleVertices = generateCircleVertices(circleRadius, numVertices);

    //// Generate ellipse vertices
    std::vector<Point> ellipseVertices = generateEllipseVertices(circleRadius, axisRatio, numVertices);

    //// Write the circle vertices to a file
    writeToFile(directoryName + "/circle_vertices.txt", circleVertices);

    //// Write the ellipse vertices to a file
    writeToFile(directoryName + "/ellipse_vertices.txt", ellipseVertices);

    //// Compute the constant-speed parametrization for the circle
    std::vector<double> params_circle = computeParametrization(circleVertices);

    //// Write the parameter values for each vertex of the circle to a file
    writeParamsToFile(directoryName + "/circle_params.txt", params_circle);

    //// Compute the constant-speed parametrization for the ellipse
    std::vector<double> params_ellipse = computeParametrization(ellipseVertices);

    //// Write the parameter values for each vertex of the ellipse to a file
    writeParamsToFile(directoryName + "/ellipse_params.txt", params_ellipse);

    //// Create the SFML window
    sf::RenderWindow window(sf::VideoMode(800, 600), "SFML Shape Drawer");

    //// Create and set up the view
    sf::View view = window.getDefaultView();
    view.setSize(static_cast<float>(window.getSize().x), static_cast<float>(window.getSize().y));

    //// Apply zoom factor
    float zoomFactor = 1.0f / 7.0f; // Zoom in by a factor of 3
    view.zoom(zoomFactor); // Note: zoom factor < 1.0 zooms in

    //// Set the view to the window
    window.setView(view);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear(sf::Color::White);

        // Draw the circle
        drawShape(window, circleVertices, sf::Color::Blue);

        // Draw the ellipse
        drawShape(window, ellipseVertices, sf::Color::Red);

        window.display();
    }

    //// Save the screenshots
    sf::Texture texture;
    texture.create(window.getSize().x, window.getSize().y);
    texture.update(window);
    sf::Image screenshot = texture.copyToImage();
    screenshot.saveToFile(directoryName + "/circle_and_ellipse.png");
}