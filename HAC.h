#ifndef HEC
#define HEC

#include "datapoints.h"
#include "CFTree.h"
#include "voronoi.h"

#include <vector>
#include <utility>
#include <queue>

double distance(const std::vector<double>& A, const std::vector<double>& B);
struct CFTree;
struct Cluster;
struct Voronoi;


//Structure for defining the dendogram
//We have history - what is happening when
//We keep the initial clusters given, but we can construct every step. This will be updated to be history dependent as 
//keeping all pointers for each active cluster is wasteful

//Keeps track of the active clusters for the dendograms. These are used there for speedup
struct ActiveClusters{
    std::vector<size_t> initial_ids;
    std::vector<double> centroid;
    size_t size;
    const size_t id;
    bool active;

    ActiveClusters(const std::vector<size_t>& initial_ids, const std::vector<double>& centroid, size_t size, size_t id):
    id(id), initial_ids(initial_ids),
    centroid(centroid), size(size), active(true){}

    //For initial once
    ActiveClusters(size_t initial, const std::vector<double>& centroid, size_t size, size_t id): 
    centroid(centroid), size(size), id(id), initial_ids({initial}), active(true){}        
};


struct ActiveClustersPQ : ActiveClusters {

    struct PQEntry {
        double dist;
        size_t other;
        bool operator>(const PQEntry& o) const { 
            return dist > o.dist; 
        }
    };

    using MinPQ = std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>>;

    MinPQ pq;
    using ActiveClusters::ActiveClusters;
};


struct Dendogram {
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

    const std::vector<const Cluster*> initial_clusters;
    std::vector<Step> history;
    size_t next_id;


    Dendogram(const std::vector<const Cluster*>& initial_clusters):
    initial_clusters(initial_clusters), next_id(initial_clusters.size()){}
    

    Cluster BuildFromNode(const ActiveClusters& act_cl) const;
    void PrintHistory(const std::string& name) const;
    void PrintSummary(const std::string& name) const;

    //Placceholder. Calling them on a this structure is meaningless. 
    //Technically simple algos can be implemented to be finding the distance every time, but its much sloweer(at least 20 times)

    virtual std::pair<size_t, size_t> FindClosest(const size_t max_threads){ return {0, 0};}

    virtual bool MergeClosest(size_t a, size_t b, const size_t max_threads){ return false;}

    //Possible closest neighbor algos
    struct PQ;

};



struct Dendogram::PQ : Dendogram {

    

    std::vector<ActiveClustersPQ> actives;

    PQ(const std::vector<const Cluster*>& initial_clusters):
    Dendogram(initial_clusters){

        const size_t n = initial_clusters.size();
        this->actives.reserve(n);

        for (size_t i = 0; i < n; i += 1) {
            this->actives.emplace_back(i, initial_clusters[i]->Centroid(), initial_clusters[i]->points.size(), i);
        }

        //O(n**2) definition of distances. Distance is O(att)
        for (size_t i = 0; i < n; i += 1) {
            for (size_t j = 0; j < n; j += 1) {
                if (i == j) {
                    continue;
                }

                double dist = distance(this->actives[i].centroid, this->actives[j].centroid);
                this->actives[i].pq.push({dist, j});
            }
        }
    }

    std::pair<size_t, size_t> FindClosest(const size_t max_threads);
    bool MergeClosest(size_t a, size_t b, const size_t max_threads);

};


#endif