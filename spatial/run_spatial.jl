#=
Recommended command line to run:
>  julia -L SpatialEvolution.jl run_spatial.jl configs/example1
=#
export spatial_result, print_spatial_result, run_trial, writeheader, writerow
#=
type spatial_result_type
  N::Int64   # meta-population size
  num_subpops::Int64   # number of subpopulations
  num_env_subpops::Int64   # number of "subpopulations" to use for env variation.  0 means use num_subpops.
  ne::Int64  # number of emmigrants in horizontal transfer
  num_attributes::Int64  # number of attributes of a variant
  mu::Float64     # innovation rate
  ngens::Int64  # number of generations after burn-in
  horiz_select::Bool       # Whether to use selection during horzontal transfer
  circular_variation::Bool    # Whether to vary ideal values in a circular fashion
  extreme_variation::Bool    # Whether to vary ideal values by randomly choosing between high and low values
  normal_stddev::Float64  # standard deviation of normal distribution of mutation perturbations
  fitness_mean::Float64
  fitness_variance::Float64
  attribute_variance::Float64
end
=#
#include("types.jl")
  
function spatial_result( N::Int64, num_subpops::Int64, num_env_subpops::Int64, ne::Int64, num_attributes::Int64, mu::Float64, ngens::Int64, burn_in::Float64,
    horiz_select::Bool, circular_variation::Bool, extreme_variation::Bool, normal_stddev::Float64 )
  return spatial_result_type( N, num_subpops, num_env_subpops, ne, num_attributes, mu, ngens, burn_in,
    horiz_select, circular_variation, extreme_variation, normal_stddev, 0.0, 0.0, 0.0 )
end

function print_spatial_result( sr::spatial_result_type )
  println("N: ", sr.N)
  println("num_subpops: ", sr.num_subpops)
  println("num_env_subpops: ", sr.num_env_subpops)
  println("ne: ", sr.ne)
  println("num_attributes: ", sr.num_attributes)
  println("mu: ", sr.mu)
  println("normal_stddev: ", sr.normal_stddev)
  println("ngens: ", sr.ngens)
  println("burn_in: ", sr.burn_in)
  println("horiz_select: ", sr.horiz_select)
  println("circular_variation: ",sr.circular_variation)
  println("extreme_variation: ",sr.extreme_variation)
  println("fitness_mean: ", sr.fitness_mean)
  println("fitness_variance: ", sr.fitness_variance)
  println("attiribute_variance: ", sr.attribute_variance)
end
#=  Moved to run.jl
function run_trials( simname::AbstractString ) 
  stream = open("$(simname).csv","w")
  println("stream: ",stream)
  println("stddev: ",normal_stddev())
  sr = SpatialEvolution.spatial_result(N,num_subpops_list[1],num_env_subpops_list[1],ne_list[1],num_attributes, mu, ngens, horiz_select, circular_variation, extreme_variation_list[1], normal_stddev() )
  writeheader( stream, num_subpops_list, sr )
  trial = 1
  for num_subpops in num_subpops_list
    for ne in ne_list
      for extreme_variation in extreme_variation_list
        sr = SpatialEvolution.spatial_result(N,num_subpops,num_env_subpops,ne,num_attributes, mu, ngens, horiz_select, circular_variation, extreme_variation, normal_stddev() )
        run_trial(sr)
        writerow(stream,trial,sr)
        print_spatial_result( sr )
        trial += 1
      end
    end
  end
end
=#

function writeheader( stream::IO, num_subpops_list::Vector{Int64}, sr::spatial_result_type )
  param_strings = [
    "# $(string(Dates.today()))",
    "# N=$(sr.N)",
    "# num_subpops_list=$(num_subpops_list)",
    #"# num_attributes=$(sr.num_attributes)",
    "# mu=$(sr.mu)",
    "# horiz_select=$(sr.horiz_select)",
    #"# num_emmigrants=$(sr.ne)",
    "# ngens=$(sr.ngens)",
    "# circular_variation=$(sr.circular_variation)",
    #"# extreme_variation=$(sr.extreme_variation)",
    "# burn_in=$(sr.burn_in)",
    "# normal_stddev=$(sr.normal_stddev)"]
  write(stream,join(param_strings,"\n"),"\n")
  heads = [
    "num_subpops",
    "subpop_size",
    "num_emigrants",
    "ext_variation",
    "num_env_subpops",
    "num_attributes",
    "mean_fitness",
    "variance_fitness",
    "attribute_variance"]
  write(stream,join(heads,","),"\n")
end
    
function writerow( stream::IO, trial::Int64, sr::spatial_result_type )
  line = Any[sr.num_subpops,
          Int(ceil(sr.N/sr.num_subpops)),
          sr.ne,
          sr.extreme_variation,
          sr.num_env_subpops,
          sr.num_attributes,
          sr.fitness_mean,
          sr.fitness_variance,
          sr.attribute_variance]
  write(stream,join(line,","),"\n")
end

#=
if length(ARGS) == 0
  simname = "configs/example2"
else
  simname = ARGS[1]
end
include("$(simname).jl")
println("simname: ",simname)
println("simtype: ",simtype)
stream = open("$(simname).csv","w")
#println("stream: ",stream)
=#


if isdefined(:simtype)
#  (avg_means,avg_variances, avg_attr_variances) = spatial_simulation( N, m, mu, copy_err_prob, ngens, burn_in, ne, num_attributes, 
    #cdfe, idfe, 
#    vtbl, subpop_properties )
  run_trials()
  #println("avg_means: ",avg_means)
  #println("avg_variances: ",avg_variances)
  #println("avg_attr_variances: ",avg_attr_variances)
  #=
  for pl in pop_list
    println(pl)
  end
  =#
end
