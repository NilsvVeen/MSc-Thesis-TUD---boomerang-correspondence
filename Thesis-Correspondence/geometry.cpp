// geometry.cpp
#include "geometry.h"
#include <cmath>
#include <stdexcept>

// Define PI if not already defined
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
