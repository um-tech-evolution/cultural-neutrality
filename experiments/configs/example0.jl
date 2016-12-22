# Configuration for running Slatkin test on the results of multiple evolutions
@everywhere const simtype = 0    # Neutral
@everywhere const n_list     = [100 ]     # Population size list
const N_mu_list= [1.0, 2.0]   # Mutation rate as a multiple of 1.0/N
const ngens = 1       # Generations after burn-in
const popsize_multiplier = 1.0
const burn_in= 2.0    # generations of burn_in as a multiple of N
const T     = 10      # Number of trials
const slat_reps=100000 # Number of reps in ewens_montecarlo
const sig_level = 0.05   # one-sided
const t2error_sides = 1

