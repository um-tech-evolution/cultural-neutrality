# Configuration for running of multiple evolutions of power conformist
const simtype = 2    # Acerbi conformist
const N_list     = [100, 500]     # Population size list
const N_mu_list= [1.0, 10.0]   # Mutation rate as a multiple of 1.0/N
const ngens = 1       # Generations after burn-in
const burn_in= 1.0    # generations of burn_in as a multiple of N
const T     = 100      # Number of trials
const slat_reps=100000 # Number of reps in ewens_montecarlo
const acer_C_list = [-0.05, 0.0, 0.05]    # acer_C_list, 0.0 is neutral, negative is anti-conformity
const acer_topsize=10  # Size of conformity toplist
