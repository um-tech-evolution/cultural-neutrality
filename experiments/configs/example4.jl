# Configuration for running Slatkin and Watterson tests on the results of multiple evolutions
@everywhere const simtype = 2    # Acerbi conformist
const acer_C_list = [0.0, -0.05]    # acer_C_list, 0.0 is neutral, negative is anti-conformity
const acer_topsize=10  # Size of conformity toplist
const popsize_multiplier = 2.0   # popsize N is this multiple of sample size n
@everywhere const n_list     = [25, 50]     # Population size list
const acer_flag_list = [true]
const N_mu_list= [1.0, 10.0]   # Mutation rate as a multiple of 1.0/N
const cprob_list = [0.2]           # probability of conformist
const acprob_list = [0.1]          # probability of anti-conformist
const bottom_list = [true]
const topsize_list=[1]  # Size of conformity toplist
const bottomsize_list=[1]  # Size of anti-conformity bottomlist
const ngens = 1       # Generations after burn-in
const burn_in= 1.0    # generations of burn_in as a multiple of N
const T     = 100      # Number of trials
const slat_reps=100000 # Number of reps in ewens_montecarlo
const sig_level = 0.025  # two-sided
const t2error_sides = 2

