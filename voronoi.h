#ifndef VORONOI
#define VORONOI

#include "datapoints.h"
#include "CFTree.h"
#include <utility>
#include <vector>
#include <mutex>

struct CFTree;

struct Voronoi{

    //These are all cells
    struct Cell{
        std::vector<double> splitting_p;
        std::vector<const Cluster*> clusters;

        //Log for adding elems in parallel.
        std::vector<const Cluster*> add_log;
        std::mutex log_mutex;

        Cell(std::vector<double> splitting_p){
            this->splitting_p = splitting_p;
        }

        //Move cell
        Cell(Cell&& other) noexcept {
            this->splitting_p = std::move(other.splitting_p);
            this->clusters = std::move(other.clusters);
            this->add_log = std::move(other.add_log);
        }

        void AddToLog(const Cluster& cluster);
        void CommitLog();
    };

    double d; //overlap distance
    size_t c; //Number of cells
    std::vector<Cell> cells;

    Voronoi(double d, CFTree* cft) {
        this->d = d;
        this->cells.reserve(cft->Nodes.size());

        for (size_t i = 0; i < cft->Nodes.size(); i += 1) {
            this->cells.emplace_back(cft->Nodes[i].Centroid());
        }

        this->c = this->cells.size();
    }

    std::vector<const std::vector<double>*> AllSplittingPoints() const;

    void SplitClusterToLog(const Cluster& cluster);
    void SplitClustersThreads(const std::vector<Cluster>& clusters, size_t num_threads);

    void PrintSummary() const;

    //Parallel On The spliting of all Clusters
    void CommitAllLogs();

    //Running HAC on a Voronoi cell



};

#endif