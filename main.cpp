#include "datapoints.h"
#include "CFTree.h"
#include "voronoi.h"
#include "HAC.h"

#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <chrono>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <csv-file-path> [num_threads]" << std::endl;
        std::cerr << "Example path: Clustering-Datasets/01.\\ UCI/banknote.csv" << std::endl;
        return 1;
    }

    std::string filename = argv[1];

    double maxR = 0.1;

    //Make sure there is at least 1 thread
    size_t num_threads = std::max<size_t>(1, std::stoul(argv[2]));

    //By testing on the datasets, it turns out that the overhead computations for the 
    //addition of extra threads is hurting more then helping
    //It still can be used needed, but for now it is excluded
    auto dataloading_start = std::chrono::high_resolution_clock::now();
    Data data(filename);
    auto dataloading_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dataloading_time = dataloading_end - dataloading_start;
    /*
    auto dataloading_start = std::chrono::high_resolution_clock::now();
    Data data(filename, num_threads);
    auto dataloading_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dataloading_time = dataloading_end - dataloading_start;
    */


    if (data.points.empty()) {
        std::cerr << "No data loaded from: " << filename << std::endl;
        return 1;
    }

    data.PrintSummary();

    auto pts_to_clusters_start = std::chrono::high_resolution_clock::now();
    std::vector<Cluster> clusters = PointsToClusters(data);
    auto pts_to_clusters_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> pts_to_clusters_time = pts_to_clusters_end - pts_to_clusters_start;

    std::cout << "Maximum radius for the CFTree: " << maxR << std::endl;
    std::cout << "Threads for Voronoi split: " << num_threads << std::endl;

    auto cftree_start = std::chrono::high_resolution_clock::now();
    CFTree tree(maxR);
    tree.IncludeClusters(clusters);
    auto cftree_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> cftree_time = cftree_end - cftree_start;
    tree.PrintSummary();

    double d = maxR;

    auto voronoi_ctor_start = std::chrono::high_resolution_clock::now();
    Voronoi v(d, &tree);
    auto voronoi_ctor_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> voronoi_ctor_time = voronoi_ctor_end - voronoi_ctor_start;


    auto voronoi_start = std::chrono::high_resolution_clock::now();
    v.SplitClustersThreads(clusters, num_threads);
    auto voronoi_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> voronoi_time = voronoi_end - voronoi_start;

    //v.PrintSummary();

    auto vdendo_ctor_start = std::chrono::high_resolution_clock::now();
    VoronoiDendogram voronoi_dendo(v);
    auto vdendo_ctor_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> vdendo_ctor_time = vdendo_ctor_end - vdendo_ctor_start;

    auto run_until_d_start = std::chrono::high_resolution_clock::now();
    voronoi_dendo.RunUntilD(num_threads);
    auto run_until_d_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> run_until_d_time = run_until_d_end - run_until_d_start;

    for (size_t i = 0; i < v.cells.size(); i += 1) {
        size_t active_after = 0;
        for (const auto& a : voronoi_dendo.cell_dendos[i].actives) {
            if (a.active) {
                active_after += 1;
            }
        }

        std::cout << "Cell " << i
                  << " | clusters before HAC: " << v.cells[i].clusters.size()
                  << " | clusters after HAC: " << active_after
                  << " | merges: " << voronoi_dendo.cell_dendos[i].history.size()
                  << std::endl;

        //local_dendogram.PrintHistory("Cell " + std::to_string(i));
    }

    auto collect_start = std::chrono::high_resolution_clock::now();
    std::vector<Cluster> all_clusters = voronoi_dendo.GetAllClusters(num_threads);
    auto collect_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> collect_time = collect_end - collect_start;

    std::cout << "Total local clusters collected: " << all_clusters.size() << std::endl;

    if (all_clusters.empty()) {
        std::cerr << "No clusters collected after Voronoi phase." << std::endl;
        return 1;
    }

    std::vector<const Cluster*> all_cluster_ptrs;

    for (size_t i = 0; i < all_clusters.size(); i += 1) {
        all_cluster_ptrs.push_back(&all_clusters[i]);
    }

    auto global_start = std::chrono::high_resolution_clock::now();
    Dendogram::PQ global_dendogram(all_cluster_ptrs);
    while (global_dendogram.MergeClosest(num_threads)) {}

    auto global_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> global_time = global_end - global_start;

    std::cout << "Global dendogram after Voronoi:" << std::endl;
    global_dendogram.PrintSummary("Global");
    //global_dendogram.PrintHistory("Global");

    //std::cout << "It took " << dataloading_time.count() << "secs to standartize the data!" << std::endl;
    std::cout << "\n=== Timing summary ===" << std::endl;
    std::cout << "Data loading:               " << dataloading_time.count() << " s" << std::endl;
    std::cout << "PointsToClusters:           " << pts_to_clusters_time.count() << " s" << std::endl;
    std::cout << "CFTree build:               " << cftree_time.count() << " s" << std::endl;
    std::cout << "Voronoi ctor:               " << voronoi_ctor_time.count() << " s" << std::endl;
    std::cout << "It took " << voronoi_time.count() << "secs to finish splitting the clusters!" << std::endl;
    std::cout << "VoronoiDendogram ctor:      " << vdendo_ctor_time.count() << " s" << std::endl;
    std::cout << "RunUntilD (parallel HAC):   " << run_until_d_time.count() << " s" << std::endl;
    std::cout << "GetAllClusters (parallel):  " << collect_time.count() << " s" << std::endl;
    std::cout << "Global HAC:                 " << global_time.count() << " s" << std::endl;

    return 0;
}