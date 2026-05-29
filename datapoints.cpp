#include "datapoints.h"

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>



// Get all fields split from a single datapoint
std::vector<std::string> split(const std::string& line, char delim) {
    std::vector<std::string> parts;
    std::string item;
    std::stringstream ss(line);

    while (std::getline(ss, item, delim)) {
        parts.push_back(item);
    }

    return parts;
}

//Implementation of the reading of the data file and putting 
//it in Data - the class corresponding to datapoint

void Data::Initialize(const std::string& filename){
    //Try to open the file, return if an issue
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Could not open file: " << filename << std::endl;
        return;
    }

    char delim = ',';

    //Get the classes' names
    std::string header;
    std::getline(file, header);

    std::vector<std::string> column_names = split(header, delim);
    
    const int n = column_names.size();

    this->column_names = std::move(column_names);

    //Get all data
    std::string line;
    int rowNumber = 1;

    std::vector<Point> points;
    //For getting coords
    std::vector<double> maxes(n, -INFINITY);
    std::vector<double> mins(n, INFINITY);

    while (std::getline(file, line)) {
        //If a typo line, continue
        if (line.empty()) {
            continue;
        }

        //Getting the value and turning them to doubles
        std::vector<std::string> values = split(line, delim);

        //we have assumed that the last value is the label
        int label_idx = n - 1;

        std::vector<double> numbers;
        int label;

        for (int i = 0; i < n; i++) {
            double val = std::stod(values[i]);
            if (i == label_idx) {
                label = val;
            } else {
                numbers.push_back(val);
            }
            if(maxes[i] < val){
                maxes[i] = val;
            }else if(mins[i] > val){
                mins[i] = val;
            }
        }
        Point p(numbers, label);
        points.push_back(p);
        rowNumber++;
    }
    this->points = std::move(points);
    this->num_classes = n - 1;
    this->max_val = std::move(maxes);
    this->min_val = std::move(mins);

}
