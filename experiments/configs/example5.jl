# Configuration for running Watterson test on the results of multiple evolutions
# run simulatios with popsize 1.0/fract_size, then sample fract_size of population
@everywhere const simtype = 1    # Power conformity
const popsize_multiplier = 2.0  # popsize N is this multiple of sample size n
@everywhere const n_list  = [50, 200]     # Population size list
const N_mu_list= [1.0, 8.0]   # Mutation rate as a multiple of 1.0/N
const ngens = 1       # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
const T     = 100      # Number of trials
const cp_list = [0.0, 0.08 ]
const cprob_list = [0.2]           # probability of conformist
const acprob_list = [0.0]          # probability of anti-conformist
const cpower_list = [0.0, 0.05]    # cpower values.  0.0 is neutral.  0.05 is positive
const acpower_list = [-0.05]       # anti-conformist power
const sig_level = 0.05   # one-sided
const slat_reps=100000 # Number of reps in ewens_montecarlo
const t2error_sides = 1
