# Example Configuration for nearly neutral infinite sites simuation
# Command line to run:  julia n_neutral.jl configs/in_example1
# Creates output file configs/in_example1.csv

const nn_simtype = 2    # nearly neutral infinite alleles
type_str = "N"
mu_list_flag =true
@everywhere const N_list     = [100,200]   # sample size, popsize N=  popsize_multiplier*n
const mu_list= [0.01,0.02]   # Mutation rate 
const L = 1           # number of loci # not used in infinite alleles
const ngens = 1000      # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
dfe=dfe_neutral
dfe_str = "dfe_neutral"

