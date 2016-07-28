# Configuration for running Watterson test on the results of multiple evolutions
const simtype = 1    # Power conformity
const N_list     = [100, 500]     # Population size list
const N_mu_list= [1.0, 10.0]   # Mutation rate as a multiple of 1.0/N
const ngens = 1       # Generations after burn-in
const burn_in= 1.0    # generations of burn_in as a multiple of N
const T     = 100      # Number of trials
const cp_list = [0.0, 0.05]
