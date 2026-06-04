#ifndef HEC_OPTIMIZE
#define HEC_OPTIMIZE

#include "datapoints.h"
#include "CFTree.h"
#include "voronoi.h"
#include "HAC.h"

#include <vector>
#include <utility>
#include <queue>

struct VoronoiDendogramLocal{ 
    
    const Voronoi& voro;
    std::vector<Dendogram::PQ> cell_dendos;
    const size_t max_threads;

    VoronoiDendogramLocal(const Voronoi& voro, const size_t max_threads):
    max_threads(max_threads), 
    voro(voro){
        this->cell_dendos.reserve(voro.cells.size());

        for (size_t i = 0; i < voro.cells.size(); i += 1) {
            this->cell_dendos.emplace_back(voro.cells[i].clusters);
        }
    }


    void RunUntilD();
    std::vector<Cluster> GetAllClusters();

    void PrintVoronoiSummary();
    
};


//Paper idea

struct VoronoiDendogramPaper{
    
     struct MergeStep {
        size_t cell;       
        size_t left_id;
        size_t right_id;   
        double dist;
    };
    std::vector<MergeStep> history;
    const Voronoi& voro;
    std::vector<Dendogram::PQ> cell_dendos;
    
    VoronoiDendogramPaper(const Voronoi& voro): 
    voro(voro){
        this->cell_dendos.reserve(voro.cells.size());

        for (size_t i = 0; i < voro.cells.size(); i += 1) {
            this->cell_dendos.emplace_back(voro.cells[i].clusters);
        }
    }
    struct OverallRes{
        double dist; 
        size_t cell; 
        size_t elem_a; 
        size_t elem_b;
    };
   
    OverallRes FindClosestOverall(const size_t max_threads);

    void DeactivateEverywhere(const Cluster* p, size_t skip_cell);
    bool Merge(const OverallRes& r);
    void RunUntilD(size_t max_threads);
    std::vector<Cluster> GetAllClustersPaper(size_t max_threads);
    

};

#endif 