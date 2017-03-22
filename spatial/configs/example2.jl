# Configuration for running spatial simulation
#=
Recommended command line to run:
>  julia -L SpatialEvolution.jl run.jl configs/example2
or
>  julia -p 4 -L SpatialEvolution.jl run.jl configs/example2
=#
export simtype
#global simtype=1 # means spatial structure by fitness or selection coefficient adjustment
# simtype=2 means spatial structure by changing the ideal values for attributes
@everywhere simtype = 2    
const T     = 40      # Number of trials
const vtbl = Dict{Int64,variant_type}()
@everywhere const N = 8        # Meta-population size
const num_subpops_list = [1,2]
const normal_stddev = 0.1
const mu = 0.05                 # per-individual innovation rate 
#const ne = 0                    # number emmigrants
const ne_list = [0,2]                    # number emmigrants
const num_attributes = 2        # number attributes for quantitative representation
const ngens = 5       # Generations after burn-in
#const use_fit_locations=false
const horiz_select=true
#const horiz_select_list=[false,true]
const circular_variation=true
#const extreme_variation=false
const extreme_variation_list=[false,true]
const burn_in= 0.1    # generations of burn_in as a multiple of N
const use_fit_locations_list=[false,true]
const ideal_max = 0.8
const ideal_min = 0.2
const ideal_range = 0.1


