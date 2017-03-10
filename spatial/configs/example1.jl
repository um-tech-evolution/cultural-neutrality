# Configuration for running spatial simulation
idfe() =  1.0+rand(Distributions.Gamma(innov_alpha,innov_beta))
cdfe() =  -rand(Distributions.Gamma(copy_alpha,copy_beta))
@everywhere const simtype = 1    # First spatial simulation type
const T     = 40      # Number of trials
const vtbl = Dict{Int64,variant_type}()
@everywhere const N = 20        # Meta-population size
const m = 4                     # Number of subpopulations
const mu = 0.05                 # per-individual innovation rate
const copy_err_prob = 0.1       # per-individual copy error probability
const ne = 1                    # number emmigrants
const ngens = 20       # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
const innov_alpha=1.0
const innov_beta=0.02
const copy_alpha=0.2
const copy_beta=0.1
