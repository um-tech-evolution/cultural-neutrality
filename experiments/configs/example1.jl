# Configuration for running Slatkin and Watterson tests on the results of multiple evolutions
const simtype = 1    # Power conformity
const T     = 3      # Number of trials
const popsize_multiplier = 2.0  # popsize N is this multiple of sample size n
@everywhere const n_list     = [25]     # sample size n list, popsize N may be larger
const N_mu_list= [2.0]   # Mutation rate as a multiple of 1.0/N
const ngens = 1       # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
cprob_list = [0.2]           # probability of conformist
acprob_list = [0.0]          # probability of anti-conformist
const cpower_list = [0.0, 0.05]    # cpower values.  0.0 is neutral.  0.05 is positive
const acpower_list = [-0.05]       # anti-conformist power
const slat_reps=100000 # Number of reps in ewens_montecarlo
const sig_level = 0.05   # one-sided
const t2error_sides = 1
