# Example Configuration for nearly neutral simuation
# Command line to run:  julia n_neutral.jl configs/nn_example
# Creates output file configs/nn_example.csv

const nn_simtype = 1    # nearly neutral infinite alleles
@everywhere const N_list     = [40,160]     # popsize 
const N_mu_list= [1.0,2.0]   # Mutation rate as a multiple of 1.0/N
const L = 3           # number of loci # not used in infinite alleles
const ngens = 1000      # Generations after burn-in
const fix_minimum = 0.6  # Fraction of population needed to "fix"
const burn_in= 2.0    # generations of burn_in as a multiple of N
dfe=dfe_neutral
dfe_str = "dfe_neutral"

