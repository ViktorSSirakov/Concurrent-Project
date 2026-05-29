#include "datapoints.h"
#include "CFTree.h"
#include "voronoi.h"

#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    //Bad Call of the executable
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <csv-file-path>" << std::endl;
        std::cerr << "Example path: Clustering-Datasets/01.\\ UCI/banknote.csv" << std::endl;        
        return 1;
    }

    //Getting the filename
    std::string filename = argv[1];

    //Initializing the data
    Data data;
    data.Initialize(filename);

    //No data was actually saved
    if (data.points.empty()) {
        std::cerr << "No data loaded from: " << filename << "\n";
        return 1;
    }


    
    //Testing couts
    std::cout << "Loaded file: " << filename << "\n";
    std::cout << "Columns: " << data.column_names.size() << "\n";
    std::cout << "Points: " << data.points.size() << "\n";
    std::cout << "Features per point: " << data.points[0].data.size() << "\n";
    for(size_t i = 0; i < data.max_val.size(); i += 1){
        std::cout << "My max values is " << data.max_val[i] << ", and the minimum is " 
        << data.min_val[i] << " for the attribute " << data.column_names[i] << std::endl;
    }
    std::cout << "\n\n\n" << std::endl;
    //Testing the CFstruct
    double maxR;
    if(argc == 3){
        std::cout << "Maximum radious for the CFTree is provided and it is " 
        << argv[2] << std::endl;

        maxR = std::stod(argv[2]);
    }else{
        maxR = 1;
    }
    CFTree tree(maxR);
    tree.IncludeData(data);

    std::vector<std::vector<double>> centroids = tree.GetAllCentroid();

    std::cout << "Points: " << data.points.size() << std::endl;
    std::cout << "maxR: " << maxR << std::endl;
    std::cout << "CF nodes / centroids: " << centroids.size() << std::endl;

    for (size_t i = 0; i < centroids.size(); i++) {
        std::cout << "Centroid " << i << ": ";

        for (size_t j = 0; j < centroids[i].size(); j++) {
            std::cout << centroids[i][j];

            if (j + 1 < centroids[i].size()) {
                std::cout << ", ";
            }
        }

        std::cout << "\n";
    }

    //Testing Voronoi 
    double d = 0.3;
    Voronoi v(d, &tree);
    v.SplitData(data);
    for(size_t i = 0; i < v.cells.size(); i += 1){
        std::cout << "Cell number " << (i + 1) << " has " << v.cells[i].points.size() << " inside." << std::endl;
    }
    return 0;
}