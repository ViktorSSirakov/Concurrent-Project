#include "HAC.h"
#include "CFTree.h"
#include "datapoints.h"
#include "voronoi.h"

#include <vector>
#include <cmath>
#include <utility>
#include <iostream>
#include <thread>
#include <atomic>


double distance(const std::vector<double>& A, const std::vector<double>& B);
struct Cluster;
// Distance between 2 clusters



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



void Dendogram::PrintHistory(const std::string& name) const {
    std::cout << "\n" << name << " history:" << std::endl;

    for (size_t i = 0; i < this->history.size(); i += 1) {
        const Dendogram::Step& step = this->history[i];

        std::cout << "  Step " << i
                  << ": " << step.left_id
                  << " + " << step.right_id
                  << " -> " << step.new_id
                  << " at distance " << step.distance
                  << std::endl;
    }
}

void Dendogram::PrintSummary(const std::string& name) const {
    std::cout << name << " merges: " << this->history.size() << std::endl;
    std::cout << name << " active clusters: " << this->active_ids.size() << std::endl;
}


static void RunCellsUntilDHelp(VoronoiDendogram* self, std::atomic<size_t>* next, double d) {

    const size_t total = self->cell_dendos.size();
    for (;;) {
        size_t i = next->fetch_add(1);
        if (i >= total) return;

        Dendogram& dendo = self->cell_dendos[i];
        while (dendo.active_ids.size() > 1) {
            std::pair<size_t, size_t> closest = dendo.FindClosest();
            double dist = distance(dendo.nodes[closest.first].centroid, dendo.nodes[closest.second].centroid);
            if (dist > d){ 
                break;
            }
            dendo.MergeClosest();
        }
    }
}

void VoronoiDendogram::RunUntilD(const size_t max_threads) {

    const size_t n = this->voro.c;
    if (n == 0) return;

    const size_t num_threads = std::max<size_t>(1, std::min(max_threads, n));

    const double d = this->voro.d;

    std::atomic<size_t> next(0);
    std::vector<std::thread> threads;

    for (size_t t = 0; t < num_threads; t += 1) {
        threads.emplace_back(RunCellsUntilDHelp, this, &next, d);
    }

    for (auto& th : threads){ 
        th.join();
    }
}

//Building the remaining clusters. I arrange allocate space for the cluster and then a 
//thread is tasked with construction. When done, it make a different one. There is an atomic variable 
//to check which cluster were working at

void CollectClustersWorker(VoronoiDendogram* self, std::atomic<size_t>* next, 
    const std::vector<size_t>* offsets, std::vector<Cluster>* out) {

    const size_t n_cells = self->cell_dendos.size();
    for (;;) {
        size_t i = next->fetch_add(1, std::memory_order_relaxed);
        if (i >= n_cells) return;

        const Dendogram& dendo = self->cell_dendos[i];
        const size_t base = (*offsets)[i];
        for (size_t j = 0; j < dendo.active_ids.size(); j += 1) {
            (*out)[base + j] = dendo.BuildFromNode(dendo.active_ids[j]);
        }
    }
}

std::vector<Cluster> VoronoiDendogram::GetAllClusters(const size_t max_threads) {
    const size_t n_cells = this->cell_dendos.size();
    if (n_cells == 0){ 
        std::cerr << "There was problem with loading" << std::endl;
        return {};
    }

    std::vector<size_t> offsets(n_cells);
    size_t total = 0;

    for (size_t i = 0; i < n_cells; i += 1) {
        offsets[i] = total;
        total += this->cell_dendos[i].active_ids.size();
    }

    std::vector<Cluster> out(total);
    if (total == 0) return out;

    const size_t num_threads = std::max<size_t>(1, std::min(max_threads, n_cells));
    std::atomic<size_t> next(0);

    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    for (size_t t = 0; t < num_threads; t += 1) {
        threads.emplace_back(CollectClustersWorker, this, &next, &offsets, &out);
    }

    for (auto& th : threads){ 
        th.join();
    }

    return out;
}