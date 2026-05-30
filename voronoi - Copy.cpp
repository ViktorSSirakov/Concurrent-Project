#include "voronoi.h"
#include "HAC.h"
#include "CFTree.h"
#include <cmath>

double distance(const std::vector<double>& A, const std::vector<double>& B);

//All spliting points we have. Not useful
std::vector<const std::vector<double>*> Voronoi::AllSplittingPoints() const{
    std::vector<const std::vector<double>*> cs;

    for (const Cell& cell : this->cells) {
        cs.push_back(&cell.splitting_p);
    }

    return cs;
}


//Split a single cluster
void Voronoi::SplitCluster(const Cluster& cluster){

    std::vector<double> centroid = cluster.Centroid();

    double min_dist = INFINITY;
    //Find best distance
    for(size_t i = 0; i < this->cells.size(); i += 1){
        double dist_cell = distance(centroid, this->cells[i].splitting_p);

        if(dist_cell < min_dist){
            min_dist = dist_cell;
        }
    }

    //Put in all cells where it should be
    for(size_t i = 0; i < this->cells.size(); i += 1){
        if(distance(centroid, this->cells[i].splitting_p) <= min_dist + this->d){
            this->cells[i].clusters.push_back(&cluster);
        }
    }

}


//Split all clusters
void Voronoi::SplitClusters(const std::vector<Cluster>& clusters){
    for(const Cluster& cluster : clusters){
        this->SplitCluster(cluster);
    }
}

void Voronoi::PrintSummary() const {
    for (size_t i = 0; i < this->cells.size(); i += 1) {
        std::cout << "Cell number " << (i + 1)
                  << " has " << cells[i].clusters.size()
                  << " inside." << std::endl;
    }
}