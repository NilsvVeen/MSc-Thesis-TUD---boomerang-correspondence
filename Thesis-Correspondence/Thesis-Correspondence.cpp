// Thesis-Correspondence.cpp : Defines the entry point for the application.
//

#include "Thesis-Correspondence.h"

// main.cpp
#include "geometry.h"
#include "file_utils.h"
#include "sfml_utils.h"

int main()
{
    // Define parameters using integers
    int circleRadius = 10;       // Radius of the circle in units (adjusted for text size)
    int numVertices = 100;       // Number of vertices on the circle and ellipse
    int axisRatio = 2;           // Ratio a/b for the ellipse (integer value)

    // Create a directory for output files
    std::string directoryName = "output_files";
    createDirectory(directoryName);

    // Generate circle vertices
    std::vector<Point> circleVertices = generateCircleVertices(circleRadius, numVertices);

    // Generate ellipse vertices
    std::vector<Point> ellipseVertices = generateEllipseVertices(circleRadius, axisRatio, numVertices);

    // Write the circle vertices to a file
    writeToFile(directoryName + "/circle_vertices.txt", circleVertices);

    // Write the ellipse vertices to a file
    writeToFile(directoryName + "/ellipse_vertices.txt", ellipseVertices);

    // Compute the constant-speed parametrization for the circle
    std::vector<double> params_circle = computeParametrization(circleVertices);

    // Write the parameter values for each vertex of the circle to a file
    writeParamsToFile(directoryName + "/circle_params.txt", params_circle);

    // Compute the constant-speed parametrization for the ellipse
    std::vector<double> params_ellipse = computeParametrization(ellipseVertices);

    // Write the parameter values for each vertex of the ellipse to a file
    writeParamsToFile(directoryName + "/ellipse_params.txt", params_ellipse);

    // Create the SFML window
    sf::RenderWindow window(sf::VideoMode(800, 600), "SFML Shape Drawer");

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

    // Save the screenshots
    sf::Texture texture;
    texture.create(window.getSize().x, window.getSize().y);
    texture.update(window);
    sf::Image screenshot = texture.copyToImage();
    screenshot.saveToFile("circle_and_ellipse.png");

    return 0;
}
