# Configuration for running Slatkin test on the results of multiple evolutions
const simtype = 1    # Slatkin test
const N_list     = [100, 500]     # Population size list
const mu_list= [1.0, 10.0]   # Mutation rate as a multiple of 1.0/N
const ngens = 1       # Generations after burn-in
const burn_in= 1.0    # generations of burn_in as a multiple of 1/mu
const T     = 10      # Number of trials
const slat_reps=100000 # Number of reps in ewens_montecarlo
