Commands to run the Algorithms:

Template 
./main <Path-To-DataSet> <num_thr> <min_datapoints> <app>

num_thr - number of threads to be used
min_datapoints - at how many cluster to continue to Phase 2 
Path_To_DataSet - The folder Clustering-Dataset is provided with example datasets. Only .csv files can be processed
app - the approach to be used. 0 is for normal HAC, 1 and 2 match the approach


Example command to run the experiment on the dataset banknote with 8 threads and the approach described in the pPOP paper - the settings used for debugging.

./main "Clustering-Datasets/01. UCI/banknote.csv" 8 10000 1


Given Commands to run the entire experiment:

datasets=(
  "Clustering-Datasets/01. UCI/banknote.csv"
  "Clustering-Datasets/01. UCI/biodeg.csv"
  "Clustering-Datasets/01. UCI/EEG Eye State.csv"
)

thr_ct=(1 8)

mkdir -p build &&
for ds in "${datasets[@]}"; do
  name=$(basename "$ds" .csv)        # e.g. "banknote"
  for t in "${thr_ct[@]}"; do        # thread counts
    for mode in 0 1 2; do
      for i in {1..5}; do
        ./main "$ds" "$t" 10000 "$mode"
      done > "build/out_${name}_t${t}_m${mode}.txt"
    done
  done
done


The command creates a folder build and fills it up with the results. 
By changing the dataset variable, you can run the same experiment on custom datasets.
The experiment is done for 1 and 8 threads. You can run it on different amounts by editing the variable thr_ct.



