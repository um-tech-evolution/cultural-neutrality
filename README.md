# Neutral Cultural Evolution Simulation

This repository contains Julia code for a paper on tests for cultural neutrality.
We intend to implement several tests for neutrality that can be applied to cultural evolution.
We will determine the sensitivity of the tests to deviations from the assumptions of the
Wright-Fisher infinite alleles model.  As of this writing, June 6, 2016, this is
"work in progress".

To run an example of the Slatkin monte carlo neutrality test appled to multiple 
runs of the Wright-Fisher infinite alleles model, using the following command line from 
the "experiments" subdirectory:

[wright@evotech1 experiments]$ julia -p 4 -L simulation.jl run.jl configs/example

For this run, the parameter settings are in the file configs/example.jl.  The results
of the run are in the CSV file example.csv.

These results can be analyzed/summarized by running the analysis.jl program, also in the
experiments subdirectory as with the following command line:

[wright@evotech1 experiments]$ julia analysis.jl configs/example
spdf: 4x6 DataFrames.DataFrame
| Row | N_mu | N   | mean_prob | median_prob | q05 | q95 |
|-----|------|-----|-----------|-------------|-----|-----|
| 1   | 1.0  | 100 | 0.49824   | 0.452585    | 0   | 10  |
| 2   | 1.0  | 500 | 0.418947  | 0.458045    | 0   | 10  |
| 3   | 10.0 | 100 | 0.639068  | 0.727705    | 0   | 10  |
| 4   | 10.0 | 500 | 0.490435  | 0.518105    | 0   | 10  |
thetadf slatkin: 4x5 DataFrames.DataFrame
| Row | N_mu | N   | mean_theta | median_theta | true_theta |
|-----|------|-----|------------|--------------|------------|
| 1   | 1.0  | 100 | 2.20546    | 1.86664      | 2.0        |
| 2   | 1.0  | 500 | 2.27483    | 2.20467      | 2.0        |
| 3   | 10.0 | 100 | 23.7804    | 23.1207      | 20.0       |
| 4   | 10.0 | 500 | 22.4472    | 22.8473      | 20.0       |


The program writes out the two shown data frames to stdout and to two CSV files
example_slat_prob.csv and example_slat_theta.csv.

By changing the simtype constant in the configuration file, the program runs the 
Watterson test instead of the Slatkin test, or can estimate K values.

We plan to add a power-law test and a turnover test.
