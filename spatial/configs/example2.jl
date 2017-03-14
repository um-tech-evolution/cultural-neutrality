# Configuration for running spatial simulation
export simtype
idfe() =  1.0+rand(Distributions.Gamma(innov_alpha,innov_beta))
cdfe() =  -rand(Distributions.Gamma(copy_alpha,copy_beta))
#global simtype=1 # means spatial structure by fitness or selection coefficient adjustment
# simtype=2 means spatial structure by changing the ideal values for attributes
@everywhere simtype = 2    
const T     = 40      # Number of trials
const vtbl = Dict{Int64,variant_type}()
@everywhere const N = 8        # Meta-population size
const num_subpops_list = [2,4]
const mu = 0.05                 # per-individual innovation rate, not used quantitative
const copy_err_prob = 1.0       # per-individual copy error probability, not used quantitative
const ne = 0                    # number emmigrants
const num_attributes = 2        # number attributes for quantitative representation
const ngens = 5       # Generations after burn-in
const burn_in= 0.1    # generations of burn_in as a multiple of N
normal_stddev() = return 0.1
const innov_alpha=1.0
const innov_beta=0.02
const copy_alpha=0.2
const copy_beta=0.1
