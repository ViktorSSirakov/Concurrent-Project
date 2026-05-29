#ifndef HEC
#define HEC

#include "datapoints.h"
#include <vector>

struct Cluster {

    std::vector<const Point*> points;
    std::vector<int> point_ids;    
    
    Cluster(const Point& p){
        this->points.push_back(&p);
        this->point_ids.push_back(p.id);
    }

    std::vector<double> Centroid() const;
    void Merge(const Cluster& other);

};
//Turning the points into clusters
std::vector<Cluster> PointsToClusters(const Data& data);

//Merging clusters with centroids distances M <= d
std::vector<Cluster> HACClustersUntilD(const std::vector<const Cluster*>& cluster_ptrs, double d);
#endif