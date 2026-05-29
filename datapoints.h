#ifndef POINTSH
#define POINTSH

#include <vector>
#include <string>


//All unique points 
struct Point{
    std::vector<double> data;
    double label;
    
    Point(std::vector<double> data, double label){
        this->data = data;
        this->label = label;
    }

};

//Structure covering all retrieved data
struct Data{
    std::vector<std::string> column_names;
    std::vector<Point> points;
    int num_classes = 0;
    std::vector<double> min_val;
    std::vector<double> max_val;

    void Initialize(const std::string& filename);
};

#endif