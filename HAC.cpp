#include "HAC.h"
#include <vector>
#include "CFTree.h"
#include <cmath>
double distance(const std::vector<double>& A, const std::vector<double>& B);

// Distance between 2 clusters
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

void Cluster::Merge(const Cluster& other){
    this->points.insert(this->points.end(), other.points.begin(), other.points.end());
    point_ids.insert(point_ids.end(), other.point_ids.begin(), other.point_ids.end());
}

//Turning all Points into clustors
std::vector<Cluster> PointsToClusters(const Data& data){
    std::vector<Cluster> clusters;

    for(const Point& p : data.points){
        clusters.push_back(Cluster(p));
    }

    return clusters;
}

std::vector<Cluster> HACClustersUntilD(const std::vector<const Cluster*>& cluster_ptrs, double d){
    std::vector<Cluster> clusters;

    //Copy the clusters from the cell so we can merge them locally
    for(const Cluster* c : cluster_ptrs){
        clusters.push_back(*c);
    }

    //Run HAC until no pair is closer than d
    while(clusters.size() > 1){
        double best_dist = INFINITY;
        size_t best_i = 0;
        size_t best_j = 1;

        //Find the closest pair of clusters
        for(size_t i = 0; i < clusters.size(); i += 1){
            for(size_t j = i + 1; j < clusters.size(); j += 1){
                double dist = ClusterDistance(clusters[i], clusters[j]);

                if(dist < best_dist){
                    best_dist = dist;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        //Stop local HAC when the closest pair is not inside distance d
        if(best_dist >= d){
            break;
        }

        //Merge closest pair
        clusters[best_i].Merge(clusters[best_j]);

        //Remove the merged-away cluster
        clusters.erase(clusters.begin() + best_j);
    }

    return clusters;
}