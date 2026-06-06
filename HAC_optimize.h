#ifndef HEC_OPTIMIZE
#define HEC_OPTIMIZE

#include "datapoints.h"
#include "CFTree.h"
#include "voronoi.h"
#include "HAC.h"

#include <vector>
#include <utility>
#include <queue>

//Attempt to fix the initialization. The function should work no matter the method used. It devides all vals into the correct place/ cell
//It outputs the dendos for the cells.
std::vector<Dendogram::PQ> BuildCellDendosParallel(const Voronoi& voro, size_t max_threads);

struct VoronoiDendogramLocal{ 
    
    const Voronoi& voro;
    std::vector<Dendogram::PQ> cell_dendos;
    const size_t max_threads;

    VoronoiDendogramLocal(const Voronoi& voro, const size_t max_threads):
    max_threads(max_threads), 
    voro(voro){
        //fixing the initialization issue
        this->cell_dendos = BuildCellDendosParallel(voro, max_threads);
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
    const size_t max_threads;

    VoronoiDendogramPaper(const Voronoi& voro, const size_t max_threads):
    max_threads(max_threads), 
    voro(voro){
        //Fixing the initialization issue
        this->cell_dendos = BuildCellDendosParallel(voro, max_threads);
    }
    
    struct OverallRes{
        double dist; 
        size_t cell; 
        size_t elem_a; 
        size_t elem_b;
    };
   
    OverallRes FindClosestOverall();

    void DeactivateEverywhere(const Cluster* p, size_t skip_cell);
    bool Merge(const OverallRes& r);
    void RunUntilD();
    std::vector<Cluster> GetAllClusters();
    

};

#endif 