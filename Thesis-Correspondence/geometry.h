// geometry.h
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <string>

struct Point {
    double x, y;
};

double distance(const Point& a, const Point& b);
std::vector<double> computeParametrization(const std::vector<Point>& vertices);
std::vector<Point> generateCircleVertices(int radius, int numVertices);
std::vector<Point> generateEllipseVertices(int circleRadius, int axisRatio, int numVertices);

void old_1(std::string& directoryName);

#endif // GEOMETRY_H
