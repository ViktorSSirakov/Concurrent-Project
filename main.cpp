#include "datapoints.h"

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

    return 0;
}