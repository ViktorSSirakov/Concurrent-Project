#ifndef HEC
#define HEC

#include "datapoints.h"
#include <vector>
#include <utility>


struct Cluster {
    std::vector<const Point*> points;
    std::vector<int> point_ids;  
    
    Cluster(const Point& p){
        this->points.push_back(&p);
        this->point_ids.push_back(p.id);
    }

    std::vector<double> Centroid() const;

};
Cluster Merge(const Cluster& first, const Cluster& second);

//Turning the points into clusters
std::vector<Cluster> PointsToClusters(const Data& data);

//Structure for defining the dendogram
//We have history - what is happening when
//We keep the initial clusters given, but we can construct every step
struct Dendogram {
    struct Node{
        std::vector<size_t> initial_ids;
        std::vector<double> centroid;
        size_t size;

        Node(const std::vector<size_t>& initial_ids, const std::vector<double>& centroid, size_t size){
            this->initial_ids = initial_ids;
            this->centroid = centroid;
            this->size = size;
        }
        //For initial once
        Node(size_t initial, const std::vector<double>& centroid, size_t size){
            this->initial_ids.push_back(initial);
            this->centroid = centroid;
            this->size = size;
        }

    };

    struct Step {
        size_t left_id;
        size_t right_id;
        size_t new_id;
        double distance;

        Step(size_t left_id, size_t right_id, size_t new_id, double distance) {
            this->left_id = left_id;
            this->right_id = right_id;
            this->new_id = new_id;
            this->distance = distance;
        }
    };

    std::vector<const Cluster*> initial_clusters;
    std::vector<Node> nodes;
    std::vector<size_t> active_ids;
    std::vector<Step> history;
    size_t next_id;

    Dendogram(const std::vector<const Cluster*>& initial_clusters) {
        size_t n = initial_clusters.size();

        for(size_t i = 0; i < n; i += 1) {
            this->initial_clusters.push_back(initial_clusters[i]);
            Node node(i, initial_clusters[i]->Centroid(), initial_clusters[i]->points.size());
            this->nodes.push_back(node);
            this->active_ids.push_back(i);
        }

        this->next_id = n;
    }

    std::pair<size_t, size_t> FindClosest() const;
    bool MergeClosest();
    Cluster BuildFromNode(size_t node_id) const;
};

struct VoronoiDendogram : Dendogram {
    double d;

    VoronoiDendogram(const std::vector<const Cluster*>& clusters, double d):
          Dendogram(clusters),
          d(d) {}

    void RunUntilD();
};

#endif