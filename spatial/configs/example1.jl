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
@everywhere const N = 8        # Meta-population size
const num_subpops_list = [2]                     # Number of subpopulations
const mu = 0.05                 # per-individual innovation rate 
#const ne = 1                    # number emmigrants
const ne_list=[1]
const num_attributes = 1        # number attributes for quantitative representation
#const use_fit_locations=true
const use_fit_locations_list=[false,true]
const horiz_select=false
const circular_variation=false
#const extreme_variation=true
const extreme_variation_list=[true]
const ngens = 2       # Generations after burn-in
const burn_in= 0.1    # generations of burn_in as a multiple of N
const normal_stddev = 0.05
const ideal_max = 0.8
const ideal_min = 0.2
const ideal_range = 0.1
