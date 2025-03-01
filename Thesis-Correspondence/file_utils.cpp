// file_utils.cpp
#include "file_utils.h"
#include <fstream>
#include <iostream>
#include <sys/stat.h> // For mkdir (POSIX) or _mkdir (Windows)
#include <direct.h>  // For _mkdir on Windows
#include <filesystem> // For clearing the folder (C++17)

#include <filesystem>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "stl_utils.h"

// Function to create a directory if it does not exist
void createDirectory(const std::string& dirName) {
#ifdef _WIN32
    _mkdir(dirName.c_str());
#else
    mkdir(dirName.c_str(), 0755);
#endif
}

// Function to clear the folder of all files
void clearDirectory(const std::string& dirName) {
    try {
        // Check if the directory exists
        if (std::filesystem::exists(dirName)) {
            // Iterate over all files in the directory and remove them
            for (const auto& entry : std::filesystem::directory_iterator(dirName)) {
                std::filesystem::remove_all(entry.path());
            }
        }
    }
    catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error clearing directory: " << e.what() << std::endl;
    }
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


#include <regex>
#include <iostream>
#include <filesystem>
#include <vector>
#include <Eigen/Dense>

void loadMeshesFromFolder(const std::string& folderPath,
    std::vector<Eigen::MatrixXd>& inputV,
    std::vector<Eigen::MatrixXi>& inputF) {

    std::regex shapePattern(R"(^Shape.*\.obj$)");  // Matches "ShapeXYZ.obj", "Shape1.obj", etc.

    try {
        for (const auto& entry : std::filesystem::directory_iterator(folderPath)) {
            std::string filename = entry.path().filename().string();

            // Debugging: Print filenames
            std::cout << "Checking file: " << filename << std::endl;

            // Check if the file is an .obj and matches "Shape*.obj"
            if (entry.path().extension() == ".obj" && std::regex_match(filename, shapePattern)) {
                std::cout << "Loading mesh: " << filename << std::endl;

                Eigen::MatrixXd V;
                Eigen::MatrixXi F;

                // Check if reading the file is successful
                if (readMeshFromFile(entry.path().string(), V, F)) {
                    inputV.push_back(V);
                    inputF.push_back(F);
                }
                else {
                    std::cerr << "Error reading mesh: " << filename << std::endl;
                }
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    }
}


void loadMeshes(const std::vector<std::string>& filePaths,
    std::vector<Eigen::MatrixXd>& inputV,
    std::vector<Eigen::MatrixXi>& inputF) {
    for (const std::string& filePath : filePaths) {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        readMeshFromFile(filePath, V, F);
        inputV.push_back(V);
        inputF.push_back(F);
    }
}
