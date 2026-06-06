#include "datapoints.h"
#include "CFTree.h"
#include "voronoi.h"
#include "HAC.h"
#include "HAC_optimize.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

struct IterationTiming {
    double maxR;
    size_t clusters_before;
    size_t clusters_after;
    size_t clusters_disappeared;

    double cftree_time = 0.0;
    double voronoi_ctor_time = 0.0;
    double voronoi_time = 0.0;
    double vdendo_ctor_time = 0.0;
    double run_until_d_time = 0.0;
    double collect_time = 0.0;

    double OverallTime(){
        const double time_cost =
            this->cftree_time +
            this->voronoi_ctor_time +
            this->voronoi_time +
            this->vdendo_ctor_time +
            this->run_until_d_time +
            this->collect_time;

        return time_cost;
    }
};

template <typename Dendogram>
std::vector<Cluster> RunLocalPhase(Voronoi& v, size_t max_threads, IterationTiming& timing) {
    auto ctor_start = std::chrono::high_resolution_clock::now();
    Dendogram voronoi_dendo(v, max_threads);
    auto ctor_end = std::chrono::high_resolution_clock::now();
    timing.vdendo_ctor_time = std::chrono::duration<double>(ctor_end - ctor_start).count();

    std::cout << "the split happened!" << std::endl;

    auto run_start = std::chrono::high_resolution_clock::now();
    voronoi_dendo.RunUntilD();
    auto run_end = std::chrono::high_resolution_clock::now();
    timing.run_until_d_time = std::chrono::duration<double>(run_end - run_start).count();

    auto collect_start = std::chrono::high_resolution_clock::now();
    std::vector<Cluster> result = voronoi_dendo.GetAllClusters();
    auto collect_end = std::chrono::high_resolution_clock::now();
    timing.collect_time = std::chrono::duration<double>(collect_end - collect_start).count();

    return result;
}

int main(int argc, char* argv[]) {
    const size_t num_args = 5;


    if (argc < num_args) {
        std::cerr << "Usage: " << argv[0] << " <csv-file-path> <num_threads> <max_cells> <approach>" << std::endl;
        std::cerr << "Example path: Clustering-Datasets/01.\\ UCI/banknote.csv 4 100 <1>" << std::endl;
        return 1;
    }


    const size_t max_datapoints = 8000;
    std::string filename = argv[1];

    size_t max_threads = std::max<size_t>(1, std::stoul(argv[2]));
    size_t max_cells = std::max<size_t>(2, std::stoul(argv[3]));

    double maxR = 0.1;

    auto dataloading_start = std::chrono::high_resolution_clock::now();
    Data data(filename, max_datapoints);
    auto dataloading_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dataloading_time = dataloading_end - dataloading_start;

    if (data.points.empty()) {
        std::cerr << "No data loaded from: " << filename << std::endl;
        return 1;
    }

    //data.PrintSummary();

    auto pts_to_clusters_start = std::chrono::high_resolution_clock::now();
    std::vector<Cluster> clusters = PointsToClusters(data);
    auto pts_to_clusters_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> pts_to_clusters_time = pts_to_clusters_end - pts_to_clusters_start;

    std::cout << "Threads: " << max_threads << std::endl;
    std::cout << "Target max cells: fewer than " << max_cells << std::endl;
    std::cout << "Initial singleton clusters: " << clusters.size() << std::endl;

    std::vector<Cluster> all_clusters;
    std::vector<IterationTiming> iteration_timings;
    double d = maxR;

    const size_t approach = std::stoul(argv[4]);
    
    if(approach >= 3){
        std::cerr << "This approach is unknown!" << std::endl;
        return 1;
    }

    if(approach > 0){
        do {
            std::cout << "\nRunning with delta = " << d << std::endl;

            const size_t clusters_before = clusters.size();

            IterationTiming timing;
            timing.maxR            = maxR;
            timing.clusters_before = clusters_before;

            auto cftree_start = std::chrono::high_resolution_clock::now();
            CFTree tree(maxR, 2000);
            tree.IncludeClusters(clusters);
            auto cftree_end = std::chrono::high_resolution_clock::now();
            timing.cftree_time = std::chrono::duration<double>(cftree_end - cftree_start).count();

            auto voronoi_ctor_start = std::chrono::high_resolution_clock::now();
            Voronoi v(d, &tree);
            auto voronoi_ctor_end = std::chrono::high_resolution_clock::now();
            timing.voronoi_ctor_time = std::chrono::duration<double>(voronoi_ctor_end - voronoi_ctor_start).count();
            std::cout << "Number of Voronoi cells: " << tree.Nodes.size() << std::endl;

            auto voronoi_start = std::chrono::high_resolution_clock::now();
            v.SplitClustersThreads(clusters, max_threads);
            auto voronoi_end = std::chrono::high_resolution_clock::now();
            timing.voronoi_time = std::chrono::duration<double>(voronoi_end - voronoi_start).count();

            if (approach == 2)
                all_clusters = RunLocalPhase<VoronoiDendogramLocal>(v, max_threads, timing);
            else
                all_clusters = RunLocalPhase<VoronoiDendogramPaper>(v, max_threads, timing);

            const size_t clusters_after = all_clusters.size();
            const size_t clusters_disappeared =
                clusters_before > clusters_after ? clusters_before - clusters_after : 0;

            timing.clusters_after       = clusters_after;
            timing.clusters_disappeared = clusters_disappeared;
            iteration_timings.push_back(timing);

            std::cout << "Clusters before current step: " << clusters_before << std::endl;
            std::cout << "Total local clusters collected: " << clusters_after << std::endl;
            std::cout << "Clusters disappeared in current step: " << clusters_disappeared << std::endl;

            clusters = all_clusters;
            maxR *= 1.1;
            d = maxR;

        } while (all_clusters.size() >= max_cells);

    }else{
        all_clusters = clusters;
    }

    std::vector<const Cluster*> all_cluster_ptrs;

    for (size_t i = 0; i < all_clusters.size(); i += 1) {
        all_cluster_ptrs.push_back(&all_clusters[i]);
    }

    auto global_start = std::chrono::high_resolution_clock::now();
    Dendogram::PQ global_dendogram(all_cluster_ptrs, max_threads);
    auto global_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> global_time = global_end - global_start;

    auto global_merge_start = std::chrono::high_resolution_clock::now();
    size_t left = global_dendogram.MakeDendogram();
    auto global_merge_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> global_merge_time = global_merge_end - global_merge_start;

    std::cout << "\nGlobal dendogram after Voronoi:" << std::endl;
    global_dendogram.PrintSummary("Global");


    std::cout << "\n=== Per-run timing summary ===" << std::endl;

    for (size_t i = 0; i < iteration_timings.size(); i += 1) {
        const auto& t = iteration_timings[i];

        std::cout << "\nRun " << i + 1 << " | maxR = " << t.maxR
                  << " | clusters before = " << t.clusters_before
                  << " | clusters after = " << t.clusters_after
                  << " | disappeared = " << t.clusters_disappeared
                  << std::endl;

        std::cout << "CFTree build:               " << t.cftree_time << " s" << std::endl;
        std::cout << "Voronoi ctor:               " << t.voronoi_ctor_time << " s" << std::endl;
        std::cout << "Voronoi split:              " << t.voronoi_time << " s" << std::endl;
        std::cout << "VoronoiDendogramLocal ctor:      " << t.vdendo_ctor_time << " s" << std::endl;
        std::cout << "RunUntilD:                  " << t.run_until_d_time << " s" << std::endl;
        std::cout << "GetAllClusters:             " << t.collect_time << " s" << std::endl;
    }

    std::cout << "\n\n================== Overall Timing Summary ===================" << std::endl;

    for(size_t i = 0; i < iteration_timings.size(); i += 1){
        std::cout << "For rerun number " << (i + 1)
                  << ", it took " << iteration_timings[i].OverallTime()
                  << " secs and removed "
                  << iteration_timings[i].clusters_disappeared
                  << " clusters."
                  << std::endl;
    }

    std::cout << "For the global case, initialization took "
              <<  global_time.count()
              << " secs." << std::endl;

    std::cout << "For the global case, making the dendogram took "
              <<  global_merge_time.count()
              << " secs." << std::endl;

    std::cout << "Overall time in the general case: " <<  global_time.count() + global_merge_time.count() << std::endl;

    std::cout << "We are left with " << left << " clusters" <<std ::endl;
    


    
    //baseline_hac.PrintSummary("Baseline HAC");
    
    return 0;
}



