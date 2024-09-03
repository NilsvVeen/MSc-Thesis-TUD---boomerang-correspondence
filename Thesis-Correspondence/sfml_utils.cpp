// sfml_utils.cpp
#include "sfml_utils.h"
#include <SFML/Graphics.hpp>

// Function to draw points in a SFML window
void drawShape(sf::RenderWindow& window, const std::vector<Point>& vertices, sf::Color color) {
    sf::VertexArray shape(sf::LineStrip, vertices.size() + 1);
    for (size_t i = 0; i < vertices.size(); ++i) {
        shape[i].position = sf::Vector2f(vertices[i].x + window.getSize().x / 2,
            -vertices[i].y + window.getSize().y / 2);
        shape[i].color = color;
    }
    shape[vertices.size()].position = sf::Vector2f(vertices[0].x + window.getSize().x / 2,
        -vertices[0].y + window.getSize().y / 2);
    shape[vertices.size()].color = color;
    window.draw(shape);
}
