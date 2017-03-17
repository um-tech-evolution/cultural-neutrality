# Configuration for running spatial simulation
#=
Recommended command line to run:
>  julia -L SpatialEvolution.jl run.jl configs/example1
or
>  julia -p 4 -L SpatialEvolution.jl run.jl configs/example1
=#
export simtype
@everywhere simtype = 2    
const T     = 40      # Number of trials
@everywhere const N = 10        # Meta-population size
const num_subpops_list = [2]                     # Number of subpopulations
const mu = 0.05                 # per-individual innovation rate, not used quantitative
#const ne = 1                    # number emmigrants
const ne_list=[1]
const num_attributes = 1        # number attributes for quantitative representation
const num_env_subpops=0
const horiz_select=false
const circular_variation=false
#const extreme_variation=true
const extreme_variation_list=[false,true]
const ngens = 5       # Generations after burn-in
const burn_in= 0.1    # generations of burn_in as a multiple of N
const normal_stddev = 0.05
