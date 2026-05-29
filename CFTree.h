#ifndef CFTREE
#define CFTREE

#include "datapoints.h"
#include <vector>
#include <iostream>
#include <cmath>


double distance(const std::vector<double>& A, const std::vector<double>& B);

//AS nothing is really saved in that datastruct and was made first with points, we keep as points, not single point clusters.
//No point in changing

struct CFNode{
    int n = 0;
    std::vector<double> ls;
    double ss = 0.0;

    CFNode(const size_t dim){
        std::vector<double> ls(dim, 0.0);
        this->ls = std::move(ls);
    }


    void AddPoint(const Point& p);
    void RemovePoint(const Point& p);
    std::vector<double> Centroid() const;
    double Radius();
};

struct CFTree{
    double maxR;
    std::vector<CFNode> Nodes;


    CFTree(double maxR){
        this->maxR = maxR; 
    }

    std::vector<std::vector<double>> GetAllCentroid();
    int ClosestEntry(const Point& p) const;
    void IncludeDatapoint(const Point& p);
    void IncludeData(const Data& data);
};


#endif