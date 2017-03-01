# Configuration for nearly neutral simuation to illustrate the Hartl & Clark
#    heterozygosity ratio formula (page 250).
# Example run from command line:  julia n_neutral.jl nn_configs/nn_example1
function dfemod(x)
  dfe_mod(x,fit_inc=1.0+sel_coef,modulus=dfe_modulus)
end
const nn_simtype = 2    # nearly neutral infinite sites
const mu_list_flag = true
#const popsize_multiplier_list = [1]  # popsize N is this multiple of sample size n
@everywhere const N_list     = [50,100]     # sample size n list, popsize N may be larger
#const N_mu_list= [2.0,4.0]   # Mutation rate as a multiple of 1.0/N
const mu_list= [0.01,0.02]   # Mutation rate as a multiple of 1.0/N
const L = 1           # number of loci
const ngens = 10      # Generations after burn-in
const fix_minimum = 0.45  # Fraction of population needed to "fix"
const burn_in= 2.0    # generations of burn_in as a multiple of 1/mu
const dfe=dfemod
const dfe_modulus = 10
const sel_coef = -0.1
const dfe_str = "dfe_mod sel coef: $(sel_coef)  modulus: $(dfe_modulus) "

