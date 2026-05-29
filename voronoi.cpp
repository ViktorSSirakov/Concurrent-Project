#include "voronoi.h"
#include "datapoints.h"
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


void Voronoi::SplitDatapoint(const Point& p){

    double min_dist = INFINITY;
    //Find best distance
    for(size_t i = 0; i < this->cells.size(); i += 1){
        double dist_cell = distance(p.data, this->cells[i].splitting_p);
        if(dist_cell < min_dist){
            min_dist = dist_cell;
        }
    }
    //Put in all cells where it should be
    for(size_t i = 0; i < this->cells.size(); i += 1){
        if(distance(p.data, this->cells[i].splitting_p) <= min_dist + this->d){
            this->cells[i].points.push_back(&p);
        }
    }

}

void Voronoi::SplitData(const Data& data){
    for(const Point& p: data.points){
        this->SplitDatapoint(p);
    }
}
