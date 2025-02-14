#include <SFML/Graphics.hpp>
#include <Eigen/Dense>
#include <vector>
#include <cmath>

// Symmetrize the curve along the Y-axis
std::vector<Eigen::Vector2d> symmetrizeCurve(const std::vector<Eigen::Vector2d>& points) {
    std::vector<Eigen::Vector2d> symmetrized;
    for (const auto& point : points) {
        symmetrized.push_back(point);                         // Original point
        symmetrized.push_back(Eigen::Vector2d(-point.x(), point.y())); // Mirrored point
    }
    return symmetrized;
}

// Smooth the curve using Chaikin's algorithm
std::vector<Eigen::Vector2d> smoothCurve(const std::vector<Eigen::Vector2d>& points, int iterations) {
    std::vector<Eigen::Vector2d> smoothed = points;
    for (int iter = 0; iter < iterations; ++iter) {
        std::vector<Eigen::Vector2d> newSmoothed;

        for (size_t i = 0; i < smoothed.size(); ++i) {
            Eigen::Vector2d p1 = smoothed[i];
            Eigen::Vector2d p2 = smoothed[(i + 1) % smoothed.size()]; // Wrap around for closed loop

            // Add 1/4 and 3/4 points
            newSmoothed.push_back(0.75 * p1 + 0.25 * p2);
            newSmoothed.push_back(0.25 * p1 + 0.75 * p2);
        }
        smoothed = newSmoothed;
    }
    return smoothed;
}

// Snap a point to the nearest grid position
Eigen::Vector2d snapToGrid(const Eigen::Vector2d& point, double gridSize) {
    return Eigen::Vector2d(
        std::round(point.x() / gridSize) * gridSize,
        std::round(point.y() / gridSize) * gridSize
    );
}

int showAndDraw() {
    sf::RenderWindow window(sf::VideoMode(800, 600), "Boomerang Drawing Tool");
    window.setFramerateLimit(60);

    std::vector<Eigen::Vector2d> drawnPoints;
    std::vector<Eigen::Vector2d> smoothedPoints;
    int smoothnessLevel = 1;
    bool isDrawing = false;

    // Reset button
    sf::RectangleShape resetButton(sf::Vector2f(150, 50));
    resetButton.setFillColor(sf::Color(200, 200, 200));
    resetButton.setOutlineThickness(2);
    resetButton.setOutlineColor(sf::Color::Black);
    resetButton.setPosition(620, 50);

    double gridSize = 10.0; // Grid size for snapping

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();

            // Handle mouse button press
            if (event.type == sf::Event::MouseButtonPressed) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    sf::Vector2i mousePos(event.mouseButton.x, event.mouseButton.y);

                    // Check if reset button is clicked
                    if (resetButton.getGlobalBounds().contains(mousePos.x, mousePos.y)) {
                        drawnPoints.clear();
                        smoothedPoints.clear();
                    }
                    else {
                        isDrawing = true;
                    }
                }
            }

            // Stop drawing on mouse release
            if (event.type == sf::Event::MouseButtonReleased) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    isDrawing = false;

                    // Smooth and symmetrize the curve
                    smoothedPoints = smoothCurve(drawnPoints, smoothnessLevel);
                    smoothedPoints = symmetrizeCurve(smoothedPoints);
                }
            }
        }

        // Capture mouse position while drawing
        if (isDrawing) {
            sf::Vector2i mousePos = sf::Mouse::getPosition(window);
            Eigen::Vector2d snappedPoint = snapToGrid(Eigen::Vector2d(mousePos.x, mousePos.y), gridSize);
            if (drawnPoints.empty() || (drawnPoints.back() - snappedPoint).norm() > gridSize * 0.5) {
                drawnPoints.push_back(snappedPoint);
            }
        }

        // Render
        window.clear(sf::Color::White);

        // Draw grid
        for (float x = 0; x < window.getSize().x; x += gridSize) {
            sf::Vertex line[] = {
                sf::Vertex(sf::Vector2f(x, 0), sf::Color(220, 220, 220)),
                sf::Vertex(sf::Vector2f(x, window.getSize().y), sf::Color(220, 220, 220))
            };
            window.draw(line, 2, sf::Lines);
        }
        for (float y = 0; y < window.getSize().y; y += gridSize) {
            sf::Vertex line[] = {
                sf::Vertex(sf::Vector2f(0, y), sf::Color(220, 220, 220)),
                sf::Vertex(sf::Vector2f(window.getSize().x, y), sf::Color(220, 220, 220))
            };
            window.draw(line, 2, sf::Lines);
        }

        // Draw raw points
        sf::VertexArray rawLine(sf::LineStrip, drawnPoints.size() + 1); // Include closing line
        for (size_t i = 0; i < drawnPoints.size(); ++i) {
            rawLine[i].position = sf::Vector2f(drawnPoints[i].x(), drawnPoints[i].y());
            rawLine[i].color = sf::Color::Black;
        }
        if (!drawnPoints.empty()) {
            rawLine[drawnPoints.size()].position = sf::Vector2f(drawnPoints[0].x(), drawnPoints[0].y()); // Close loop
            rawLine[drawnPoints.size()].color = sf::Color::Black;
        }
        window.draw(rawLine);

        // Draw smoothed and symmetrized points
        sf::VertexArray smoothLine(sf::LineStrip, smoothedPoints.size() + 1);
        for (size_t i = 0; i < smoothedPoints.size(); ++i) {
            smoothLine[i].position = sf::Vector2f(smoothedPoints[i].x(), smoothedPoints[i].y());
            smoothLine[i].color = sf::Color::Blue;
        }
        if (!smoothedPoints.empty()) {
            smoothLine[smoothedPoints.size()].position = sf::Vector2f(smoothedPoints[0].x(), smoothedPoints[0].y());
            smoothLine[smoothedPoints.size()].color = sf::Color::Blue;
        }
        window.draw(smoothLine);

        // Draw reset button
        window.draw(resetButton);

        window.display();
    }

    return 0;
}
