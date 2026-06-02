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


bool Dendogram::PQ::MergeClosest(const size_t max_threads) {

    auto [a, b] = FindClosest(max_threads);

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



//Voronoi Dendogram stuff now. 
void RunCellsUntilDHelp(VoronoiDendogram* self, std::atomic<size_t>* next, double d) {

    const size_t total = self->cell_dendos.size();
    for (;;) {
        size_t i = next->fetch_add(1);
        if (i >= total) return;

        Dendogram::PQ& dendo = self->cell_dendos[i];
        for (;;) {
            auto [a, b] = dendo.FindClosest(1);
            if (a == b) {
                break;
            }
            double dist = distance(dendo.actives[a].centroid, dendo.actives[b].centroid);
            
            if (dist > d) {
                break;
            }
            dendo.MergeClosest(1);
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

        const Dendogram::PQ& dendo = self->cell_dendos[i];
        const size_t base = (*offsets)[i];

        size_t j = 0;
        for (const auto& a : dendo.actives) {
            if (!a.active) {
                continue;
            }
            (*out)[base + j] = dendo.BuildFromNode(a);
            j += 1;
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

        size_t active_count = 0;
        for (const auto& a : this->cell_dendos[i].actives) {
            if (a.active) {
                active_count += 1;
            }
        }
        total += active_count;
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