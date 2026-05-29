#ifndef POINTSH
#define POINTSH

#include <vector>
#include <string>

double distance(const std::vector<double>& A, const std::vector<double>& B);


//All unique points 
struct Point{
    std::vector<double> data;
    double label;
    int id;
    Point(std::vector<double> data, double label, int id){
        this->data = data;
        this->label = label;
        this->id = id;
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