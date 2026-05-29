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

int main(int argc, char* argv[]) {
    if(argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <csv-file-path>" << std::endl;
        std::cerr << "Example path: Clustering-Datasets/01.\\ UCI/banknote.csv" << std::endl;
        return 1;
    }

    std::string filename = argv[1];

    Data data;
    data.Initialize(filename);

    if(data.points.empty()) {
        std::cerr << "No data loaded from: " << filename << std::endl;
        return 1;
    }

    std::cout << "Loaded file: " << filename << std::endl;
    std::cout << "Columns: " << data.column_names.size() << std::endl;
    std::cout << "Points: " << data.points.size() << std::endl;
    std::cout << "Features per point: " << data.points[0].data.size() << std::endl;

    for(size_t i = 0; i < data.max_val.size(); i += 1) {
        std::cout << "My max values is " << data.max_val[i]
                << ", and the minimum is " << data.min_val[i]
                << " for the attribute " << data.column_names[i]
                << std::endl;
    }

    std::cout << "\n\n\n" << std::endl;

    std::vector<Cluster> clusters = PointsToClusters(data);

    double maxR;
    if(argc == 3) {
        std::cout << "Maximum radious for the CFTree is provided and it is "
                << argv[2] << std::endl;
        maxR = std::stod(argv[2]);
    } else {
        maxR = 1;
    }

    CFTree tree(maxR);
    tree.IncludeClusters(clusters);

    std::vector<std::vector<double>> centroids = tree.GetAllCentroid();

    std::cout << "Points: " << data.points.size() << std::endl;
    std::cout << "maxR: " << maxR << std::endl;
    std::cout << "CF nodes / centroids: " << centroids.size() << std::endl;

    double d = 1;
    Voronoi v(d, &tree);
    v.SplitClusters(clusters);

    for(size_t i = 0; i < v.cells.size(); i += 1) {
        std::cout << "Cell number " << (i + 1)
                << " has " << v.cells[i].clusters.size()
                << " inside." << std::endl;
    }

    std::vector<Cluster> all_clusters;

    for(size_t i = 0; i < v.cells.size(); i += 1) {
        VoronoiDendogram local_dendogram(v.cells[i].clusters, d);
        local_dendogram.RunUntilD();

        std::cout << "Cell " << i
                << " | points: " << v.cells[i].clusters.size()
                << " | clusters after HAC: " << local_dendogram.active_ids.size()
                << " | merges: " << local_dendogram.history.size()
                << std::endl;

        PrintHistory(local_dendogram, "Cell " + std::to_string(i));

        for(size_t id : local_dendogram.active_ids) {
            Cluster cluster = local_dendogram.BuildFromNode(id);
            all_clusters.push_back(cluster);
        }
    }

    std::cout << "\nTotal local clusters collected: "
            << all_clusters.size()
            << std::endl;

    if(all_clusters.empty()) {
        std::cerr << "No clusters collected after Voronoi phase." << std::endl;
        return 1;
    }

    std::vector<const Cluster*> all_cluster_ptrs;

    for(size_t i = 0; i < all_clusters.size(); i += 1) {
        all_cluster_ptrs.push_back(&all_clusters[i]);
    }

    Dendogram global_dendogram(all_cluster_ptrs);

    while(global_dendogram.active_ids.size() > 1) {
        global_dendogram.MergeClosest();
    }

    std::cout << "\nGlobal dendogram after Voronoi:" << std::endl;
    std::cout << "Global merges: "
            << global_dendogram.history.size()
            << std::endl;

    PrintHistory(global_dendogram, "Global");

    std::cout << "Final active clusters: "
            << global_dendogram.active_ids.size()
            << std::endl;

    return 0;
}