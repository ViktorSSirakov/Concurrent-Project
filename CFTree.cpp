#include "CFTree.h"
#include "datapoints.h"
#include <vector>
#include <iostream>
#include <cmath>

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

//Adding an elem to the CFnode
void CFNode::AddPoint(const Point& p){
    this->n += 1;
    for(size_t i = 0; i < p.data.size(); i += 1){
        this->ls[i] += p.data[i];
        this->ss += p.data[i] * p.data[i];
    }
}

//Removing an elem 
void CFNode::RemovePoint(const Point& p){
    this->n -= 1;
    for(size_t i = 0; i < p.data.size(); i += 1){
        this->ls[i] -= p.data[i];
        this->ss -= p.data[i] * p.data[i];
    }
}

//Getting the centroid
std::vector<double> CFNode::Centroid() const{
    std::vector<double> c(this->ls.size(), 0.0);
    if(this->n <= 0){
        return c;
    }

    for (size_t i = 0; i < ls.size(); i += 1) {
        c[i] = ls[i] / n;
    }
    return c;
}

//Getting the radious of the data (needed to make sure the CF nodes arent too wide)
double CFNode::Radius(){
    if (n == 0) {
        return 0.0;
    }
    //Formula for the radious ...
    
    double lsNormSq = 0.0;
    for (double v : ls) {
        lsNormSq += v * v;
    }

    double variance = (ss / n) - (lsNormSq / (n * n));

    if (variance < 0.0) {
        variance = 0.0;
    }

    return std::sqrt(variance);
}

//CFTree from now on

std::vector<std::vector<double>> CFTree::GetAllCentroid(){
    std::vector<std::vector<double>> cs;
    for(CFNode& node: this->Nodes){
        cs.push_back(node.Centroid());
    }
    return cs;
}

//FInd closest Node.
int CFTree::ClosestEntry(const Point& p) const{
    if(this->Nodes.size() == 0){
        return -1;
    }
    double min_dist = INFINITY;
    size_t idx = 0;
    for(size_t i = 0; i < this->Nodes.size(); i += 1){
        double dist_to_point = distance(p.data, this->Nodes[i].Centroid());
        if(dist_to_point < min_dist){
            min_dist = dist_to_point;
            idx = i;
        }
    }
    return idx;
}

//Include a single point. 
//If the radious becomes too big, remove it and make a new Node
//Id nodes is empty, make a new node
void CFTree::IncludeDatapoint(const Point& p){
    int i = this->ClosestEntry(p);
    if(i == -1){
        //Handle the empty case
        CFNode cf(p.data.size());
        cf.AddPoint(p);
        this->Nodes.push_back(cf);
        return;
    }
    
    this->Nodes[i].AddPoint(p);
    
    //Handle when the radious becomes too big.
    if(this->Nodes[i].Radius() > this->maxR){
        this->Nodes[i].RemovePoint(p);
        CFNode cf(p.data.size());
        cf.AddPoint(p);
        this->Nodes.push_back(cf);
    }
    
}

//Include the whole data
void CFTree::IncludeData(const Data& data){
    for(const Point& p: data.points){
        this->IncludeDatapoint(p);
    }
}
