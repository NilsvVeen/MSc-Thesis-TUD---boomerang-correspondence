// sfml_utils.h
#ifndef SFML_UTILS_H
#define SFML_UTILS_H

#include <vector>
#include "geometry.h"
#include <SFML/Graphics.hpp>

void drawShape(sf::RenderWindow& window, const std::vector<Point>& vertices, sf::Color color);

#endif // SFML_UTILS_H
