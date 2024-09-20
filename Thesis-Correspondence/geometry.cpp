// geometry.cpp
#include "geometry.h"
#include <stdexcept>





#include "file_utils.h"

#include "geometry.h"
#include "file_utils.h"
#include "sfml_utils.h"
#include "stl_utils.h"


#ifndef PI
#define PI 3.14159265358979323846
#endif

// Function to compute the Euclidean distance between two points
double distance(const Point& a, const Point& b) {
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    return std::sqrt(dx * dx + dy * dy);
}

std::vector<double> computeParametrization(const std::vector<Point>& vertices) {
    int n = vertices.size();
    if (n < 2) {
        throw std::invalid_argument("At least two vertices are required to form a curve.");
    }

    std::vector<double> params(n, 0.0);
    double lengthCurve = 0.0;
    std::vector<double> segments(n - 1, 0.0);

    for (int i = 0; i < n - 1; ++i) {
        segments[i] = distance(vertices[i], vertices[i + 1]);
        lengthCurve += segments[i];
    }

    // Compute parameter values for each vertex
    for (int i = 1; i < n; ++i) {
        params[i] = params[i - 1] + (segments[i - 1] / lengthCurve);
    }

    // Ensure the last vertex is mapped to 1
    params[n - 1] = 1.0;

    return params;
}

// Function to generate vertices of a circle using integer values
std::vector<Point> generateCircleVertices(int radius, int numVertices) {
    std::vector<Point> vertices;
    for (int i = 0; i < numVertices; ++i) {
        double angle = 2 * PI * i / numVertices;
        double x = radius * cos(angle);
        double y = radius * sin(angle);
        vertices.push_back({ x, y });
    }
    return vertices;
}

// Function to generate vertices of an ellipse using integer values
std::vector<Point> generateEllipseVertices(int circleRadius, int axisRatio, int numVertices) {
    std::vector<Point> circleVertices = generateCircleVertices(circleRadius, numVertices);

    // Estimate the ellipse axes a and b
    int a = circleRadius * axisRatio;
    int b = circleRadius / axisRatio;

    // Generate the ellipse vertices by scaling the circle vertices
    std::vector<Point> ellipseVertices;
    for (const auto& vertex : circleVertices) {
        double x = vertex.x * a / circleRadius;
        double y = vertex.y * b / circleRadius;
        ellipseVertices.push_back({ x, y });
    }
    return ellipseVertices;
}


void old_1(std::string& directoryName) {
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