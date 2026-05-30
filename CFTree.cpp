#include "CFTree.h"
#include "datapoints.h"
#include "HAC.h"
#include <vector>
#include <iostream>
#include <cmath>



//Adding an elem to the CFnode
void CFNode::AddCluster(const Cluster& cluster){
    std::vector<double> centroid = cluster.Centroid();

    this->n += 1;

    for(size_t i = 0; i < centroid.size(); i += 1){
        this->ls[i] += centroid[i];
        this->ss += centroid[i] * centroid[i];
    }
}

//Removing an cluster. Used for removing high radious 
void CFNode::RemoveCluster(const Cluster& cluster){
    std::vector<double> centroid = cluster.Centroid();

    this->n -= 1;

    for(size_t i = 0; i < centroid.size(); i += 1){
        this->ls[i] -= centroid[i];
        this->ss -= centroid[i] * centroid[i];
    }
}

//Getting the centroid
std::vector<double> CFNode::Centroid() const{
    std::vector<double> c(this->ls.size(), 0.0);

    if(this->n <= 0){
        return c;
    }

    for(size_t i = 0; i < ls.size(); i += 1){
        c[i] = ls[i] / n;
    }

    return c;
}

//Getting the radious of the data (needed to make sure the CF nodes arent too wide)
double CFNode::Radius(){
    if(n == 0){
        return 0.0;
    }

    double lsNormSq = 0.0;
    for(double v : ls){
        lsNormSq += v * v;
    }

    double variance = (ss / n) - (lsNormSq / (n * n));

    if(variance < 0.0){
        variance = 0.0;
    }

    return std::sqrt(variance);
}

//CFTree from now on

std::vector<std::vector<double>> CFTree::GetAllCentroid(){
    std::vector<std::vector<double>> cs;

    for(CFNode& node : this->Nodes){
        cs.push_back(node.Centroid());
    }

    return cs;
}

//FInd closest Node.
int CFTree::ClosestEntry(const Cluster& cluster) const{
    if(this->Nodes.size() == 0){
        return -1;
    }

    std::vector<double> centroid = cluster.Centroid();

    double min_dist = INFINITY;
    size_t idx = 0;

    for(size_t i = 0; i < this->Nodes.size(); i += 1){
        double dist_to_cluster = distance(centroid, this->Nodes[i].Centroid());

        if(dist_to_cluster < min_dist){
            min_dist = dist_to_cluster;
            idx = i;
        }
    }

    return idx;
}

//Include a single point. 
//If the radious becomes too big, remove it and make a new Node
//Id nodes is empty, make a new node
void CFTree::IncludeCluster(const Cluster& cluster){
    int i = this->ClosestEntry(cluster);

    if(i == -1){
        CFNode cf(cluster.Centroid().size());
        cf.AddCluster(cluster);
        this->Nodes.push_back(cf);
        return;
    }

    this->Nodes[i].AddCluster(cluster);

    if(this->Nodes[i].Radius() > this->maxR){
        this->Nodes[i].RemoveCluster(cluster);

        CFNode cf(cluster.Centroid().size());
        cf.AddCluster(cluster);
        this->Nodes.push_back(cf);
    }
}

//Include the whole data
void CFTree::IncludeClusters(const std::vector<Cluster>& clusters){
    for(const Cluster& cluster : clusters){
        this->IncludeCluster(cluster);
    }
}


void CFTree::PrintSummary() const {
    std::cout << "maxR: " << this->maxR << std::endl;
    std::cout << "CF nodes / centroids: " << this->Nodes.size() << std::endl;
}
