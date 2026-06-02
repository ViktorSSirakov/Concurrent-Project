#include "datapoints.h"

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <thread>

//Distance between 2 vectors
double distance(const std::vector<double>& A, const std::vector<double>& B){
    if(A.size() != B.size()){
        std::cerr << "Dimensions mismatch!" << std::endl;
        return -1;
    }
    double dis = 0.0;
    for(size_t i = 0; i < A.size(); i += 1){
        dis += (B[i] - A[i]) * (B[i] - A[i]);
    }
    return sqrt(dis);
}


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
    int id = 0;
    while (std::getline(file, line)) {
        //If a typo line, continue
        if (line.empty()) {
            continue;
        }

        //Getting the value and turning them to doubles
        std::vector<std::string> values = split(line, delim);

        //that was a mistake. Not used for now.
        int label_idx = n + 1;
        int label;

        std::vector<double> numbers;
        

        for (int i = 0; i < n; i++) {
            double val = std::stod(values[i]);
            if (i == label_idx) {
                label = val;
            } else {
                numbers.push_back(val);
            }
            if(maxes[i] < val){
                maxes[i] = val;
            } 
            if(mins[i] > val){
                mins[i] = val;
            }
        }
        Point p(numbers, label, id);
        points.push_back(p);
        id += 1;
        rowNumber += 1;
        
    }
    this->points = std::move(points);
    this->num_classes = n - 1;
    this->max_val = std::move(maxes);
    this->min_val = std::move(mins);

}


//Test prints
// datapoints.cpp

void Data::PrintSummary() const {
    std::cout << "\n========================== Data Summary =========================" << std::endl;
    std::cout << "Loaded file: " << this->filename << std::endl;
    std::cout << "Columns: " << this->column_names.size() << std::endl;
    std::cout << "Points: " << this->points.size() << std::endl;

    for (size_t i = 0; i < max_val.size(); i += 1) {
        std::cout << "Max value is " << max_val[i]
                  << ", and the minimum is " << min_val[i]
                  << " for the attribute " << column_names[i]
                  << std::endl;
    }
}

//Standartizing values. For now, putting the values between 0 and 1. Done with multiple threads.
void StandartizationHelp(Data& data, const std::vector<double>& intervals, size_t begin, size_t end, size_t dim){
    for(size_t j = begin; j < end; j += 1){
        auto& p = data.points[j];
        for(size_t i = 0; i < dim; i += 1){
            p.data[i] /= intervals[i];
        }
    }
}

void Data::Standartize(const size_t max_threads){
    if(this->points.size() <= 0){
        std::cerr << "Empty dataset. No initialilzation was made!" << std::endl;
        return;
    }
    const size_t dim = this->points[0].data.size();
    std::vector<double> intervals(dim, 0);

    for(size_t i = 0; i < dim; i += 1){
        intervals[i] = this->max_val[i] - this->min_val[i];
        if(intervals[i] <= 0){
            std::cerr << "Something wierd is happening with the data!" << std::endl;
            return; 
        }
    }

    const size_t n = this->points.size();
    const size_t n_threads = std::max<size_t>(1, std::min(max_threads, n));

    std::vector<std::thread> threads;

    const size_t chunk = n / n_threads;
    const size_t rem   = n % n_threads;
    size_t start = 0;

    for(size_t t = 0; t < n_threads; t += 1){
        const size_t count = chunk + (t < rem ? 1 : 0);
        threads.emplace_back(StandartizationHelp, std::ref(*this), std::cref(intervals), start, start + count, dim);
        start += count;
    }
    for(auto& th : threads) th.join();
}


//Cluster stuff
double ClusterDistance(const Cluster& a, const Cluster& b){
    return distance(a.Centroid(), b.Centroid());
}


//Centroid of a cluster
std::vector<double> Cluster::Centroid() const{
    size_t dim = this->points[0]->data.size();
    std::vector<double> c(dim, 0);

    for(const Point* p: this->points){
        for(size_t i = 0; i < dim; i += 1){
            c[i] += p->data[i];
        }
    }

    for(size_t i = 0; i < dim; i += 1){
        c[i] /= this->points.size();
    }

    return c;
}

Cluster Merge(const Cluster& first, const Cluster& second){
    Cluster cl = first;
    cl.points.insert(cl.points.end(), second.points.begin(), second.points.end());
    cl.point_ids.insert(cl.point_ids.end(), second.point_ids.begin(), second.point_ids.end());
    return cl;
}

//Turning all Points into clustors
std::vector<Cluster> PointsToClusters(const Data& data){
    std::vector<Cluster> clusters;

    for(const Point& p : data.points){
        clusters.push_back(Cluster(p));
    }

    return clusters;
}

