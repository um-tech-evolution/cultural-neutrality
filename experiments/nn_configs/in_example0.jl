# Example Configuration for nearly neutral infinite sites simuation
# Run with small values of N, ngens for testing
# Command line to run:  julia n_infsites.jl configs/in_example0
# Creates output file configs/in_example0.csv

const nn_simtype = 2    # nearly neutral infinite sites
#type_str = "N"
mu_list_flag =true
@everywhere const N_list     = [5]   # sample size, popsize N=  popsize_multiplier*n
const mu_list= [0.2]   # Mutation rate 
const L = 1           # number of loci # not used in infinite alleles
const ngens = 3      # Generations after burn-in
const burn_in= 0.0    # generations of burn_in as a multiple of N
dfe=dfe_mixed
dfe_str = "dfe_mixed"

