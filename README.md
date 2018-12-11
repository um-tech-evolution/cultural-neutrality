# Neutral Cultural Evolution Simulation

Added December 11, 2018:  We have not followed up on writing these results into papers.
I have updated files to get rid of many but not all deprecation warnings for julia v6.

The simulations of infinite alleles and infinite sites has been moved to the nearly-neutral 
replository.  The infinite sites simulation does not work here because the site_collection.jl
file is missing.

The code in slatpart.jl and ewens.jl works and has been checked.

The folder simple_inf_alleles/ has been added with code for checking a computation
of equation 3.92 of Ewens (2004) for the nearly neutral paper.


From 2016:

This repository contains Julia code for a paper on tests for cultural neutrality.
We intend to implement several tests for neutrality that can be applied to cultural evolution.
We will determine the sensitivity of the tests to deviations from the assumptions of the
Wright-Fisher infinite alleles model.  As of this writing, September 13, 2016, the implemented
"deviations" from neutrality are conformist and anti-conformist copying under two models
of conformity which I call power conformity and Acerbi conformity.  Acerbi conformity
is just another name for toplist conformity.  There is also a "nearly neutral" deviation
from neutrality which may not be working correctly.

The statistical tests implemented are the Slatkin exact test, and the Watterson homozygosity
test.  

The Slatkin test is implemented in 3 ways.  Slatkin provided C code which implements
the exact test and a sampling version of the C test. Mark Madsen has somewhat refactored the
code and George Lesica has implmented calling this code in C.  The interface is in the file
src/slatkin_C.jl.  The C code is compiled into a shared library, and this library must be in
the LD_LIBRARY_PATH on a Linux machine.  Currently, the code only runs on Linux.  I have 
implemented my own "exact" version of the test in Julia, and I am working on the sampling
version.

The Watterson homozygosity is simple to compute.  Determining the significance current
requires computing all configurations for a given n and k (part of the implementation of
the Slatkin exact test).  

I have also implemented a variation of the Watterson test called p-homozygosity.  So
far, it has not proved to be promising.

To run an example of the Slatkin monte carlo neutrality test appled to multiple 
runs of the Wright-Fisher infinite alleles model, using the following command line from 
the "experiments" subdirectory:

[wright@evotech1 experiments]$ julia -p 4 -L simulation.jl run.jl configs/example1
[wright@evotech1 experiments]$ julia hyptest.jl configs/example1

For this run, the parameter settings are in the file configs/example1.jl.  The results
of the simulation are in the CSV file configs/example1.csv, and the results of the 
analysis done by hyptest.jl are in the file  configs/example1_hyptest.jl.  
The hyptest.jl run will also show the output dataframe as text output.

Additional examples are given in the test/test_simulation.sh file. (These are working
in Dec. 2018 using julia v6 with many deprecation warnings.)

Note that the RCall Julia package needs to be able to find an R installation,
and the "poweRlaw" R package needs to be installed.
If these requirements are not met, comment out the  'include("../Rcall/poweRlaw.jl")'
line in the src/NeutralCulturalEvolution.jl file.
