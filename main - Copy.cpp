#include "datapoints.h"
#include "CFTree.h"
#include "voronoi.h"
#include "HAC.h"

#include <iostream>
#include <string>
#include <vector>

void PrintHistory(const Dendogram& dendogram, const std::string& name) {
    std::cout << "\n" << name << " history:" << std::endl;

    for(size_t i = 0; i < dendogram.history.size(); i += 1) {
        const Dendogram::Step& step = dendogram.history[i];

        std::cout << "  Step " << i
                << ": " << step.left_id
                << " + " << step.right_id
                << " -> " << step.new_id
                << " at distance " << step.distance
                << std::endl;
    }
}

// main.cpp

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <csv-file-path>" << std::endl;
        std::cerr << "Example path: Clustering-Datasets/01.\\ UCI/banknote.csv" << std::endl;
        return 1;
    }

    std::string filename = argv[1];

    Data data(filename);

    if (data.points.empty()) {
        std::cerr << "No data loaded from: " << filename << std::endl;
        return 1;
    }

    data.PrintSummary();

    std::vector<Cluster> clusters = PointsToClusters(data);

    double maxR = 1;
    if (argc == 3) {
        maxR = std::stod(argv[2]);
        std::cout << "Maximum radius for the CFTree is provided and it is " << maxR << std::endl;
    }

    CFTree tree(maxR);
    tree.IncludeClusters(clusters);
    tree.PrintSummary();

    double d = 1;
    Voronoi v(d, &tree);
    v.SplitClusters(clusters);
    v.PrintSummary();

    std::vector<Cluster> all_clusters;

    for (size_t i = 0; i < v.cells.size(); i += 1) {
        VoronoiDendogram local_dendogram(v.cells[i].clusters, d);
        local_dendogram.RunUntilD();

        std::cout << "Cell " << i
                  << " | points: " << v.cells[i].clusters.size()
                  << " | clusters after HAC: " << local_dendogram.active_ids.size()
                  << " | merges: " << local_dendogram.history.size()
                  << std::endl;

        local_dendogram.PrintHistory("Cell " + std::to_string(i));

        for (size_t id : local_dendogram.active_ids) {
            Cluster cluster = local_dendogram.BuildFromNode(id);
            all_clusters.push_back(cluster);
        }
    }

    std::cout << "\nTotal local clusters collected: "
              << all_clusters.size()
              << std::endl;

    if (all_clusters.empty()) {
        std::cerr << "No clusters collected after Voronoi phase." << std::endl;
        return 1;
    }

    std::vector<const Cluster*> all_cluster_ptrs;

    for (size_t i = 0; i < all_clusters.size(); i += 1) {
        all_cluster_ptrs.push_back(&all_clusters[i]);
    }

    Dendogram global_dendogram(all_cluster_ptrs);

    while (global_dendogram.active_ids.size() > 1) {
        global_dendogram.MergeClosest();
    }

    std::cout << "\nGlobal dendogram after Voronoi:" << std::endl;
    global_dendogram.PrintSummary("Global");
    global_dendogram.PrintHistory("Global");

    return 0;
}