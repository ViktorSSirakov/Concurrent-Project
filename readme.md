Commands to run the Algorithms:

Template 
./main <Path-To-DataSet> <num_thr> <min_datapoints>

num_thr - number of threads to be used
min_datapoints - at how many cluster to continue to Phase 2 
Path_To_DataSet - The folder Clustering-Dataset is provided with example datasets. Only .csv files can be processed


Given COmmands for the experiments:

#Experiment 1
./main Clustering-Datasets/01.\ UCI/banknote.csv 1 10000
 for i in {1..5}; do ./main Clustering-Datasets/01.\ UCI/banknote.csv 1 10000; done > out.txt