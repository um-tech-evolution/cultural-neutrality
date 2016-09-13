# Configuration for running Slatkin and Watterson tests on the results of multiple evolutions
# const simtype = 1    # Power conformity
const cp_list = [0.0, -0.05]    # cpower values.  0.0 is neutral.  -0.05 is anti-conformity
const popsize_multiplier = 2.0  # popsize is this multiple of sample size
@everywhere const n_list     = [25, 50]     # sample size list
const N_mu_list= [1.0, 10.0]   # Mutation rate as a multiple of 1.0/N
const ngens = 1       # Generations after burn-in
const burn_in= 1.0    # generations of burn_in as a multiple of N
const T     = 100      # Number of trials
const slat_reps=100000 # Number of reps in ewens_montecarlo
const sig_level = 0.05   # one-sided
const t2error_sides = 1
