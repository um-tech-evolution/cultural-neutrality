# Configuration for running Slatkin and Watterson tests on the results of multiple evolutions
p_homoz_flag = true
@everywhere const simtype = 1    
const cp_list = [0.0, -0.08]    # cpower values.  0.0 is neutral.  >0.0 is positive
const popsize_multiplier = 1.0  # popsize N is this multiple of sample size n
@everywhere const n_list     = [50]     # sample size list, popsize N may be larger
const N_mu_list= [1.0, 5.0 ]   # Mutation rate as a multiple of 1.0/N
const ngens = 1       # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
const T     = 100      # Number of trials
const slat_reps=100000 # Number of reps in ewens_montecarlo
const cprob_list = [0.2]           # probability of conformist
const acprob_list = [0.0]
const cpower_list = [0.0, 0.05]    # cpower values.  0.0 is neutral.  0.05 is positive
const acpower_list = [-0.05]       # anti-conformist power
const t2error_sides = 1
# The 3 following lines ar Used to generate significance tests and are not used in simulation
const sig_level = 0.05   # one-sided
const sig_colsyms = [:w_homoz, :p_1_4]
const sig_functs = [watterson_homozygosity,p_1_4] 
