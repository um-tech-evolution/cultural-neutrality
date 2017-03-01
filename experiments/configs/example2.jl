# Configuration for running Slatkin and Watterson tests on the results of multiple evolutions
@everywhere const simtype = 2    # Acerbi conformist
#const simtype = 2    # Acerbi conformist
const n_list     = [25]     # sample size list.  Popsize is popsize_multiplier times sample size
const popsize_multiplier = 2.0  # popsize N is this multiple of sample size n
const N_mu_list= [1.0, 2.0]   # Mutation rate as a multiple of 1.0/N
const cprob_list = [0.0, 0.2]    # probability of conformity
const acprob_list = [0.20]    # probability of anti_conformity
const topsize_list=[1]  # Size of conformity toplist
const bottomsize_list=[1]  # Size of anti-conformity bottomlist
const acer_flag_list = [true]
const bottom_list = [true]
const ngens = 1       # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
const T     = 10      # Number of trials
const slat_reps=100000 # Number of reps in ewens_montecarlo
const sig_level = 0.05   # one-sided
const t2error_sides = 1
