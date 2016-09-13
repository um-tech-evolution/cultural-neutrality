# Configuration for running Slatkin and Watterson tests on the results of multiple evolutions
 #const simtype = 3    # Nearly Neutral
const popsize_multiplier = 2.0  # popsize N is this multiple of sample size n
@everywhere const n_list     = [50]     # sample size list, popsize N may be larger
const N_mu_list= [1.0, 5.0 ]   # Mutation rate as a multiple of 1.0/N
const ngens = 1       # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
const T     = 100      # Number of trials
const slat_reps=100000 # Number of reps in ewens_montecarlo
const t2error_sides = 2
const dfe_alpha = 1.0
const dfe_theta = 1.0
dfe(N, nn_select) = dfe_advantageous( N, nn_select, alpha=dfe_alpha, theta=dfe_theta )
# The 3 following lines are Used to generate significance tests and are not used in simulation
const sig_level = 0.05   # one-sided
const sig_colsyms = [:w_homoz, :p_3_0]
