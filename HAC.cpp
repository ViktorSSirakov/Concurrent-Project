#include "HAC.h"
#include <vector>
#include <cmath>
#include <utility>


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



//All about the HAC
std::pair<size_t, size_t> Dendogram::FindClosest() const {
    double best_dist = INFINITY;
    std::pair<size_t, size_t> best_pair = {-1, -1};

    for (size_t i = 0; i < active_ids.size(); i++) {
        size_t left_id = active_ids[i];

        for (size_t j = i + 1; j < active_ids.size(); j++) {
            size_t right_id = active_ids[j];

            double dist = distance(
                nodes[left_id].centroid,
                nodes[right_id].centroid
            );

            if (dist < best_dist) {
                best_dist = dist;
                best_pair = {left_id, right_id};
            }
        }
    }

    return best_pair;
}


bool Dendogram::MergeClosest() {
    if (active_ids.size() < 2) {
        return false;
    }

    std::pair<size_t, size_t> closest = FindClosest();

    size_t left_id = closest.first;
    size_t right_id = closest.second;

    const Node& left = this->nodes[left_id];
    const Node& right = this->nodes[right_id];

    double dist = distance(left.centroid, right.centroid);

    std::vector<size_t> merged_initial_ids = left.initial_ids;
    merged_initial_ids.insert(merged_initial_ids.end(), right.initial_ids.begin(), right.initial_ids.end());

    size_t merged_size = left.size + right.size;

    std::vector<double> merged_centroid(left.centroid.size(), 0.0);

    for (size_t i = 0; i < merged_centroid.size(); i += 1) {
        merged_centroid[i] = (left.centroid[i] * left.size + right.centroid[i] * right.size)/ merged_size;
    }

    size_t new_id = this->next_id;
    this->next_id += 1;

    Node merged_node(merged_initial_ids, merged_centroid, merged_size);
    this->nodes.push_back(merged_node);

    Step step(left_id, right_id, new_id, dist);
    history.push_back(step);

    for (size_t i = 0; i < active_ids.size(); i += 1) {
        if (this->active_ids[i] == left_id || this->active_ids[i] == right_id) {
            this->active_ids.erase(this->active_ids.begin() + i);
            i -= 1;
        }
    }

    this->active_ids.push_back(new_id);

    return true;
}


Cluster Dendogram::BuildFromNode(size_t node_id) const{
    const Node& node = this->nodes[node_id];
    Cluster res = *this->initial_clusters[node.initial_ids[0]];
    for(size_t i = 1; i < node.initial_ids.size(); i += 1){
        res = Merge(res, *this->initial_clusters[node.initial_ids[i]]);
    }
    return res;
}

void VoronoiDendogram::RunUntilD() {
    while(this->active_ids.size() > 1) {
        std::pair<size_t, size_t> closest = this->FindClosest();
        double dist = distance(nodes[closest.first].centroid, nodes[closest.second].centroid);
        if(dist > d) {
            return;
        }
        MergeClosest();
    }
}