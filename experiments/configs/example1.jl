# Configuration for running Slatkin and Watterson tests on the results of multiple evolutions
# const simtype = 1    # Power conformity
const cp_list = [0.0, 0.05]    # cpower values.  0.0 is neutral.  0.05 is positive
const popsize_multiplier = 2.0  # popsize N is this multiple of sample size n
@everywhere const n_list     = [25]     # sample size n list, popsize N may be larger
const N_mu_list= [2.0]   # Mutation rate as a multiple of 1.0/N
const ngens = 1       # Generations after burn-in
const burn_in= 1.0    # generations of burn_in as a multiple of N
const T     = 10      # Number of trials
const slat_reps=100000 # Number of reps in ewens_montecarlo
const sig_level = 0.05   # one-sided
const t2error_sides = 1
