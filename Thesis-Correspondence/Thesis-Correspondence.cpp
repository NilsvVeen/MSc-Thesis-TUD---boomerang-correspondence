// Thesis-Correspondence.cpp : Defines the entry point for the application.
//

#include "Thesis-Correspondence.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm> // for std::min and std::max
#include <fstream>  // For file handling
#include <sys/stat.h> // For mkdir (POSIX) or _mkdir (Windows)
#include <direct.h> // For _mkdir on Windows


#ifndef PI
#define PI 3.14159265358979323846
#endif

// Structure to represent a point (vertex) in 2D
//struct Point {
//	double x, y;
//};
struct Point {
    int x, y;
};


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

    double lengthCurve = 0.0f;
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
        int x = static_cast<int>(radius * cos(angle));
        int y = static_cast<int>(radius * sin(angle));
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
        int x = vertex.x * a / circleRadius;
        int y = vertex.y * b / circleRadius;
        ellipseVertices.push_back({ x, y });
    }
    return ellipseVertices;
}

// Function to draw points in a text grid
void drawShape(const std::vector<Point>& vertices, int width, int height) {
    std::vector<std::vector<char>> grid(height, std::vector<char>(width, ' '));

    int centerX = width / 2;
    int centerY = height / 2;

    for (const auto& vertex : vertices) {
        int x = vertex.x + centerX;
        int y = vertex.y + centerY;

        if (x >= 0 && x < width && y >= 0 && y < height) {
            grid[y][x] = '*';
        }
    }

    // Print the grid
    for (const auto& row : grid) {
        for (const auto& cell : row) {
            std::cout << cell;
        }
        std::cout << std::endl;
    }
}

// Function to create a directory if it does not exist
void createDirectory(const std::string& dirName) {
#ifdef _WIN32
    _mkdir(dirName.c_str());
#else
    mkdir(dirName.c_str(), 0755);
#endif
}

// Function to write data to a text file
void writeToFile(const std::string& filename, const std::vector<Point>& vertices) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (const auto& vertex : vertices) {
            file << vertex.x << " " << vertex.y << std::endl;
        }
        file.close();
    }
    else {
        std::cerr << "Unable to open file " << filename << std::endl;
    }
}

void writeParamsToFile(const std::string& filename, const std::vector<double>& params) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (size_t i = 0; i < params.size(); ++i) {
            file << "p_" << i + 1 << " = " << params[i] << std::endl;
        }
        file.close();
    }
    else {
        std::cerr << "Unable to open file " << filename << std::endl;
    }
}



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

    // Output number of points to a file
    std::ofstream pointsFile(directoryName + "/points_info.txt");
    if (pointsFile.is_open()) {
        pointsFile << "Amount of points in circle: " << circleVertices.size() << std::endl;
        pointsFile << "Amount of points in ellipse: " << ellipseVertices.size() << std::endl;
        pointsFile.close();
    }
    else {
        std::cerr << "Unable to open points_info.txt" << std::endl;
    }

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

    return 0;
}
