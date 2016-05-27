# Neutral Cultural Evolution Simulation

This repository contains Julia code for a paper on tests for cultural neutrality.
We intend to implement several tests for neutrality that can be applied to cultural evolution.
We will determine the sensitivity of the tests to deviations from the assumptions of the
Wright-Fisher infinite alleles model.  As of this writing, May 27, 2016, this is
"work in progress".

To run an example of the Slatkin test appled to multiple runs of the Wright-Fisher infinite alleles model,
using the following command line from the "experiments" subdirectory:

[wright@evotech1 experiments]$ julia -p 4 -L simulation.jl run.jl configs/example

For this run, the parameter settings are in the file configs/example.jl.  The results
of the run are in the CSV file example.csv.

These results can be analyzed/summarized by running the analysis.jl program, also in the
experiments subdirectory as with the following command line:

[wright@evotech1 experiments]$ julia analysis.jl configs/example
spdf: 4×6 DataFrames.DataFrame
│ Row │ N_mu │ N   │ mean_prob │ median_prob │ q025 │ q975 │
├─────┼──────┼─────┼───────────┼─────────────┼──────┼──────┤
│ 1   │ 2.0  │ 200 │ 0.591695  │ 0.626255    │ 0    │ 0    │
│ 2   │ 2.0  │ 400 │ 0.432234  │ 0.452815    │ 1    │ 10   │
│ 3   │ 5.0  │ 200 │ 0.561486  │ 0.526335    │ 0    │ 10   │
│ 4   │ 5.0  │ 400 │ 0.509664  │ 0.56529     │ 0    │ 10   │
stdf: 4×6 DataFrames.DataFrame
│ Row │ N_mu │ N   │ mean_theta │ median_theta │ q025    │ q975    │
├─────┼──────┼─────┼────────────┼──────────────┼─────────┼─────────┤
│ 1   │ 2.0  │ 200 │ 3.93393    │ 4.09041      │ 2.19562 │ 5.34354 │
│ 2   │ 2.0  │ 400 │ 4.1668     │ 4.14607      │ 2.55531 │ 6.59857 │
│ 3   │ 5.0  │ 200 │ 9.74819    │ 9.32693      │ 6.91634 │ 15.0557 │
│ 4   │ 5.0  │ 400 │ 10.8848    │ 11.068       │ 7.90883 │ 13.1308 │

The program writes out the two shown data frames to stdout and to two CSV files
example_slat_prob.csv and example_slat_theta.csv.
