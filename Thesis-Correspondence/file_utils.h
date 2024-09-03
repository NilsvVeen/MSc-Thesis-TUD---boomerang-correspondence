// file_utils.h
#ifndef FILE_UTILS_H
#define FILE_UTILS_H

#include <vector>
#include <string>
#include "geometry.h"

void createDirectory(const std::string& dirName);
void writeToFile(const std::string& filename, const std::vector<Point>& vertices);
void writeParamsToFile(const std::string& filename, const std::vector<double>& params);

#endif // FILE_UTILS_H
