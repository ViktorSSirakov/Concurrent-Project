#ifndef VORONOI
#define VORONOI

#include "CFTree.h"
#include "datapoints.h"
#include <vector>



struct Voronoi{

    //These are all cells
    struct Cell{
        std::vector<double> splitting_p;
        std::vector<const Point*> points;
        Cell(std::vector<double> splitting_p){
            this->splitting_p = splitting_p;
        } 
    };

    double d; //overlap distance
    size_t c; //Number of cells
    std::vector<Cell> cells;

    Voronoi(double d, CFTree* cft){
        this->d = d;
        for(size_t i = 0; i < cft->Nodes.size(); i += 1){
            Cell c(cft->Nodes[i].Centroid());
            this->cells.push_back(c);
        }
        c = this->cells.size();
    }

    std::vector<const std::vector<double>*> AllSplittingPoints() const;
    void SplitDatapoint(const Point& p);
    void SplitData(const Data& data);

};

#endif