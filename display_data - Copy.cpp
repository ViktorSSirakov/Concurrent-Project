//This is a dummy file to explore the data!!!


#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

struct Point {
    double x;
    double y;
    string label;
};

vector<string> split(string line, char delim) {
    vector<string> out;
    string item;
    stringstream ss(line);
    while (getline(ss, item, delim)) {
        out.push_back(item);
    }
    return out;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input-csv-path> [output-svg-path]\n";
        return 1;
    }

    string inputFile = argv[1];
    string outputFile = argc >= 3 ? argv[2] : "uci_plot.svg";

    int xCol = 0;      // first feature
    int yCol = 1;      // second feature
    int labelCol = 4;  // last column in banknote.csv

    //If no file found, report an issue!
    ifstream in(inputFile);
    if (!in) {
        cerr << "Could not open " << inputFile << "\n";
        return 1;
    }

    string line;
    getline(in, line); // skip CSV header

    vector<Point> points;

    //GO through all points and get the first and second  colomn - the only once we actually plot here
    //For banknote, column 4 is label.
    while (getline(in, line)) {
        if (line.empty()) continue;

        char delim = line.find(';') != string::npos ? ';' : ',';
        vector<string> cols = split(line, delim);

        if ((int)cols.size() <= max({xCol, yCol, labelCol})) continue;

        Point p;
        p.x = stod(cols[xCol]);
        p.y = stod(cols[yCol]);
        p.label = cols[labelCol];
        points.push_back(p);
    }

    //No points in the file!
    if (points.empty()) {
        cerr << "No points loaded.\n";
        return 1;
    }

    //For the plotting, find the bundaries of the field.
    double minX = points[0].x, maxX = points[0].x;
    double minY = points[0].y, maxY = points[0].y;

    for (const auto& p : points) {
        minX = min(minX, p.x);
        maxX = max(maxX, p.x);
        minY = min(minY, p.y);
        maxY = max(maxY, p.y);
    }

    int width = 900;
    int height = 700;
    int pad = 60;

    vector<string> colors = {
        "#2563eb", "#dc2626", "#16a34a", "#9333ea",
        "#f59e0b", "#0891b2", "#be123c", "#4b5563"
    };

    map<string, string> labelColors;
    int colorIndex = 0;

    ofstream out(outputFile);
    out << "<svg xmlns='http://www.w3.org/2000/svg' width='" << width
        << "' height='" << height << "'>\n";
    out << "<rect width='100%' height='100%' fill='white'/>\n";

    //GIve color to the datapoint
    for (const auto& p : points) {
        if (!labelColors.count(p.label)) {
            labelColors[p.label] = colors[colorIndex % colors.size()];
            colorIndex++;
        }
        //Place each point on the plot.
        double sx = pad + (p.x - minX) / (maxX - minX) * (width - 2 * pad);
        double sy = height - pad - (p.y - minY) / (maxY - minY) * (height - 2 * pad);

        out << "<circle cx='" << sx << "' cy='" << sy
            << "' r='4' fill='" << labelColors[p.label] << "' opacity='0.75'/>\n";
    }

    out << "</svg>\n";

    cout << "Wrote " << outputFile << " with " << points.size() << " points.\n";
    return 0;
}