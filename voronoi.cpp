#include "voronoi.h"
#include "HAC.h"
#include "CFTree.h"

#include <cmath>
#include <thread>


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
void Voronoi::SplitClusterToLog(const Cluster& cluster){
    std::vector<double> centroid = cluster.Centroid();

    double min_dist = INFINITY;

    for(size_t i = 0; i < this->cells.size(); i += 1){
        double dist_cell = distance(centroid, this->cells[i].splitting_p);

        if(dist_cell < min_dist){
            min_dist = dist_cell;
        }
    }

    for(size_t i = 0; i < this->cells.size(); i += 1){
        if(distance(centroid, this->cells[i].splitting_p) <= min_dist + this->d){
            this->cells[i].AddToLog(cluster);
        }
    }
}

//This is the spliting into Voronoi
void SplitClustersThreadWorker(Voronoi* voronoi, const std::vector<Cluster>& clusters, size_t start, size_t end){

    for(size_t i = start; i < end; i += 1){
        voronoi->SplitClusterToLog(clusters[i]);
    }
}

void Voronoi::SplitClustersThreads(const std::vector<Cluster>& clusters, size_t num_threads){
    if(num_threads > clusters.size()){
        num_threads = clusters.size();
    }

    std::vector<std::thread> threads;
    size_t chunk_size = (clusters.size() + num_threads - 1) / num_threads;

    for(size_t t = 0; t < num_threads; t += 1){
        size_t start = t * chunk_size;
        size_t end = start + chunk_size;

        if(end > clusters.size()){
            end = clusters.size();
        }

        threads.push_back(std::thread(SplitClustersThreadWorker, this, std::ref(clusters), start, end));
    }

    for(std::thread& thread : threads){
        thread.join();
    }

    this->CommitAllLogs();
}




void Voronoi::PrintSummary() const {
    for (size_t i = 0; i < this->cells.size(); i += 1) {
        std::cout << "Cell number " << (i + 1)
                  << " has " << cells[i].clusters.size()
                  << " inside." << std::endl;
    }
}

//Commiting all logs for the Voronoi.
void Voronoi::CommitAllLogs(){
    for(Cell& cell : this->cells){
        cell.CommitLog();
    }
}

//Cell commiting of the log
void Voronoi::Cell::AddToLog(const Cluster& cluster){
    std::lock_guard<std::mutex> lock(this->log_mutex);
    this->add_log.push_back(&cluster);
}

void Voronoi::Cell::CommitLog(){
    this->clusters.insert(
        this->clusters.end(),
        this->add_log.begin(),
        this->add_log.end()
    );

    this->add_log.clear();
}