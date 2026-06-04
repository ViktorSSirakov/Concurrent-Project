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


Cluster Dendogram::BuildFromNode(const ActiveClusters& act_cl) const {
    Cluster result;

    for (size_t id : act_cl.initial_ids) {
        const Cluster* c = initial_clusters[id];
        for (const auto& p : c->points) {
            result.points.push_back(p);
        }
    }

    return result;
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
    std::cout << name << " initial clusters: " << this->initial_clusters.size() << std::endl;
}


//PQ functions
//For the Voronoi cells, we already use our threads for different cells, but for the remainder we can use more threads!
using thread_res = std::pair<std::pair<size_t, size_t>, double>;;
void FindClosestThreadsHelper(Dendogram::PQ* dend, size_t begin, size_t end, thread_res* out) {

    double best_dist = INFINITY;
    size_t best_a = 0;
    size_t best_b = 0;

    for (size_t i = begin; i < end; i += 1) {
        if (!dend->actives[i].active){
            continue;
        }

        auto& pq = dend->actives[i].pq;

        //Pop the inactives lazily!
        while (!pq.empty() && !dend->actives[pq.top().other].active) {
            pq.pop();
        }

        if (pq.empty()){
            continue;
        }

        if (pq.top().dist < best_dist) {
            best_dist = pq.top().dist;
            best_a = i;
            best_b = pq.top().other;
        }
    }
    *out = {{best_a, best_b}, best_dist};
}

std::pair<size_t, size_t> Dendogram::PQ::FindClosest(const size_t max_threads) {


    size_t n = this->actives.size();
    const size_t num_threads = std::max<size_t>(1, std::min(max_threads, n));
    const size_t chunk = (n + num_threads - 1) / num_threads;
    std::vector<thread_res> thr_res;
    thr_res.reserve(num_threads);
    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    for (size_t t = 0; t < num_threads; t += 1) {
        const size_t begin = t * chunk;
        const size_t end = std::min(begin + chunk, n);
        if (begin >= end) {
            break;
        }
        threads.emplace_back(&FindClosestThreadsHelper, this, begin, end, &thr_res[t]);
    }

    for (auto& th : threads) {
        th.join();
    }
    double best_dist = thr_res[0].second;
    size_t idx = 0;
    for(size_t i = 1; i < thr_res.size(); i += 1){
        const double dist = thr_res[i].second;
        if(best_dist > dist){
            best_dist = dist;
            idx = i;
        }
    }
    return thr_res[idx].first;
}


bool Dendogram::PQ::MergeClosest(size_t a, size_t b, const size_t max_threads) {


    if (a == b) {
        return false;
    }

    auto& A = actives[a];
    auto& B = actives[b];

    if (!A.active || !B.active) {
        return false;
    }

    const size_t new_size = A.size + B.size;

    std::vector<double> new_centroid(A.centroid.size(), 0.0);

    for (size_t i = 0; i < new_centroid.size(); i += 1) {
        new_centroid[i] = (A.size * A.centroid[i] + B.size * B.centroid[i]) / new_size;
    }

    std::vector<size_t> new_initial_ids = A.initial_ids;
    new_initial_ids.insert(new_initial_ids.end(), B.initial_ids.begin(), B.initial_ids.end());

    const size_t new_id = next_id++;

    history.emplace_back(A.id, B.id, new_id, distance(A.centroid, B.centroid));

    A.active = false;
    B.active = false;

    actives.emplace_back(new_initial_ids, new_centroid, new_size, new_id);

    size_t new_idx = actives.size() - 1;

    for (size_t i = 0; i < new_idx; i += 1) {

        if (!actives[i].active) {
            continue;
        }

        double d = distance(actives[i].centroid, new_centroid);

        actives[i].pq.push({d, new_idx});
        actives[new_idx].pq.push({d, i});
    }

    return true;
}



