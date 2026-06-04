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
    const std::string filename;
    const size_t max_lines;


    Data(const std::string& filename, const size_t max_lines):
    max_lines(max_lines),
    filename(filename){
        this->Initialize(filename);
        this->Standartize(1);
    }

    Data(const std::string& filename, const size_t max_threads, const size_t max_lines):
    max_lines(max_lines),
    filename(filename){
        this->Initialize(filename);
        this->Standartize(max_threads);
    }


    void Initialize(const std::string& filename);
    void Standartize(const size_t num_threads);
    void PrintSummary() const;
};

struct Cluster {
    std::vector<const Point*> points;
    std::vector<int> point_ids;
    
    Cluster() = default;
    Cluster(const Point& p){
        this->points.push_back(&p);
        this->point_ids.push_back(p.id);
    }
    std::vector<double> Centroid() const;

};
Cluster Merge(const Cluster& first, const Cluster& second);

//Turning the points into clusters
std::vector<Cluster> PointsToClusters(const Data& data);





#endif