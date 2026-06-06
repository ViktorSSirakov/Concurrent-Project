#include "HAC.h"
#include "CFTree.h"
#include "datapoints.h"
#include "voronoi.h"
#include "HAC_optimize.h"

#include <vector>
#include <cmath>
#include <utility>
#include <iostream>
#include <thread>
#include <atomic>


//Voronoi Dendogram stuff now. 
void RunCellsUntilDHelp(VoronoiDendogramLocal* self, std::atomic<size_t>* next, double d) {

    const size_t total = self->cell_dendos.size();
    for (;;) {
        size_t i = next->fetch_add(1);
        if (i >= total) return;

        Dendogram::PQ& dendo = self->cell_dendos[i];
        for (;;) {
            auto [a, b] = dendo.FindClosest();
            if (a == b) {
                break;
            }
            double dist = distance(dendo.actives[a].centroid, dendo.actives[b].centroid);
            
            if (dist > d) {
                break;
            }
            dendo.MergeClosest(a, b);
        }
    }
}

void VoronoiDendogramLocal::RunUntilD() {

    const size_t n = this->voro.c;
    if (n == 0) return;

    const size_t num_threads = std::max<size_t>(1, std::min(this->max_threads, n));

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

void CollectClustersWorker(VoronoiDendogramLocal* self, std::atomic<size_t>* next,
    const std::vector<size_t>* offsets, std::vector<Cluster>* out) {

    const size_t n_cells = self->cell_dendos.size();
    for (;;) {
        size_t i = next->fetch_add(1);
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

std::vector<Cluster> VoronoiDendogramLocal::GetAllClusters() {
    
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

    const size_t num_threads = std::max<size_t>(1, std::min(this->max_threads, n_cells));
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


void VoronoiDendogramLocal::PrintVoronoiSummary(){
    std::cout << "\n\n\n============================ Voronoi Dendogram Info ===========================" << std::endl;
    for (size_t i = 0; i < this->voro.cells.size(); i += 1) {
        size_t active_after = 0;
        for (const auto& a : this->cell_dendos[i].actives) {
            if (a.active) {
                active_after += 1;
            }
        }

        std::cout << "Cell " << i
                  << " | clusters before HAC: " << this->voro.cells[i].clusters.size()
                  << " | clusters after HAC: " << active_after
                  << " | merges: " << this->cell_dendos[i].history.size()
                  << std::endl;

    }
    std::cout << std::endl; 
}





// Voronoi Dendogram Paper


void FindClossest(VoronoiDendogramPaper* self, std::atomic<size_t>* next, VoronoiDendogramPaper::OverallRes* out) {
    const size_t total = self->cell_dendos.size();
    double best = INFINITY;
    size_t bcell = 0, ba = 0, bb = 0;

    while(true){
        const size_t i = next->fetch_add(1);
        if (i >= total) break;

        Dendogram::PQ& dendo = self->cell_dendos[i];

        for (size_t k = 0; k < dendo.actives.size(); k += 1) {
            if (!dendo.actives[k].active) continue;

            auto& pq = dendo.actives[k].pq;
            while (!pq.empty() && !dendo.actives[pq.top().other].active) pq.pop();
            if (pq.empty()) continue;

            const double dist = pq.top().dist;
            if (dist < best){ 
                best = dist; 
                bcell = i; 
                ba = k; 
                bb = pq.top().other;
            }
        }
    }
    *out = {best, bcell, ba, bb};
}

VoronoiDendogramPaper::OverallRes VoronoiDendogramPaper::FindClosestOverall() {
    const size_t n = this->cell_dendos.size();
    if (n == 0) return {INFINITY, 0, 0, 0};

    const size_t num_threads = std::max<size_t>(1, std::min(this->max_threads, n));
    std::atomic<size_t> next(0);

    std::vector<OverallRes> reses(num_threads, {INFINITY, 0, 0, 0});
    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    for (size_t t = 0; t < num_threads; t += 1) {
        threads.emplace_back(FindClossest, this, &next, &reses[t]);
    }
    for (auto& th : threads) th.join();

    OverallRes best = {INFINITY, 0, 0, 0};
    for (const auto& r : reses) {
        if(best.dist > r.dist) best = r;
    }
    return best;
}


void VoronoiDendogramPaper::DeactivateEverywhere(const Cluster* p, size_t skip_cell) {
    for (size_t i = 0; i < this->cell_dendos.size(); i += 1) {

        if (i == skip_cell) continue;

        Dendogram::PQ& dendo = this->cell_dendos[i];
        for (size_t k = 0; k < dendo.actives.size(); k += 1) {
            auto& node = dendo.actives[k];
            if (!node.active || node.initial_ids.size() != 1) continue;

            if (dendo.initial_clusters[node.initial_ids[0]] == p) {
                node.active = false;                    
                break; 
            }
        }
    }
}

bool VoronoiDendogramPaper::Merge(const OverallRes& r) {
    if (r.dist >= this->voro.d) return false;

    Dendogram::PQ& home = cell_dendos[r.cell];
    if (r.elem_a == r.elem_b || !home.actives[r.elem_a].active || !home.actives[r.elem_b].active) return false;

    const auto& na = home.actives[r.elem_a];
    if (na.initial_ids.size() == 1) DeactivateEverywhere(home.initial_clusters[na.initial_ids[0]], r.cell);

    const auto& nb = home.actives[r.elem_b];
    if (nb.initial_ids.size() == 1) DeactivateEverywhere(home.initial_clusters[nb.initial_ids[0]], r.cell);

    history.push_back({r.cell, home.actives[r.elem_a].id, home.actives[r.elem_b].id, r.dist});
    home.MergeClosest(r.elem_a, r.elem_b);
    return true;
}

void VoronoiDendogramPaper::RunUntilD() {
    while (true) {
        OverallRes r = FindClosestOverall();
        if (!Merge(r)) break;
    }
}

void CollectClustersWorkerPaper(VoronoiDendogramPaper* self, std::atomic<size_t>* next,
    const std::vector<size_t>* offsets, std::vector<Cluster>* out) {

    const size_t n_cells = self->cell_dendos.size();
    for (;;) {
        size_t i = next->fetch_add(1);
        if (i >= n_cells) return;

        const Dendogram::PQ& dendo = self->cell_dendos[i];
        const size_t base = (*offsets)[i];

        size_t j = 0;
        for (const auto& a : dendo.actives) {
            if (!a.active) continue;
            (*out)[base + j] = dendo.BuildFromNode(a);
            j += 1;
        }
    }
}

std::vector<Cluster> VoronoiDendogramPaper::GetAllClusters() {
    const size_t n_cells = this->cell_dendos.size();
    if (n_cells == 0) {
        std::cerr << "There was problem with loading" << std::endl;
        return {};
    }

    std::vector<size_t> offsets(n_cells);
    size_t total = 0;
    for (size_t i = 0; i < n_cells; i += 1) {
        offsets[i] = total;
        size_t active_count = 0;
        for (const auto& a : this->cell_dendos[i].actives)
            if (a.active) active_count += 1;
        total += active_count;
    }

    std::vector<Cluster> out(total);
    if (total == 0) return out;

    const size_t num_threads = std::max<size_t>(1, std::min(this->max_threads, n_cells));
    std::atomic<size_t> next(0);

    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (size_t t = 0; t < num_threads; t += 1) {
        threads.emplace_back(CollectClustersWorkerPaper, this, &next, &offsets, &out);
    }
    for (auto& th : threads) th.join();

    return out;
}


