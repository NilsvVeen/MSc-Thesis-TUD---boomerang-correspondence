// file_utils.cpp
#include "file_utils.h"
#include <fstream>
#include <iostream>
#include <sys/stat.h> // For mkdir (POSIX) or _mkdir (Windows)
#include <direct.h>  // For _mkdir on Windows

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
