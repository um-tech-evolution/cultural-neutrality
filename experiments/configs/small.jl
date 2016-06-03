const N_list   = [100]   # Population size list
const mu_list    = [1.0, 5.0, 10.0]   # Mutation rate as a multiple of  1.0/N
const ngens = 1       # Generations
const burn_in= 2.0    # generations of burn_in as a multiple of N
const T     = 100    # Number of trials
const slat_reps=5000 # slat_reps == 0 means use Watterson instead of Slatkin
