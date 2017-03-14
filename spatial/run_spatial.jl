#=
Recommended command line to run:
>  julia -L SpatialEvolution.jl run_spatial.jl configs/example1
=#
type spatial_result_type
  N::Int64   # meta-population size
  num_subpops::Int64   # number of subpopulations
  ne::Int64  # number of emmigrants in horizontal transfer
  num_attributes::Int64  # number of attributes of a variant
  ngens::Int64  # number of generations after burn-in
  normal_stddev::Float64  # standard deviation of normal distribution of mutation perturbations
  fitness_mean::Float64
  fitness_variance::Float64
  attribute_variance::Float64
end
  
function spatial_result( N::Int64, num_subpops::Int64, ne::Int64, num_attributes::Int64, ngens::Int64, normal_stddev::Float64 )
  return spatial_result_type( N, num_subpops, ne, num_attributes, ngens, normal_stddev, 0.0, 0.0, 0.0 )
end

function print_spatial_result( sr::spatial_result_type )
  println("N: ", sr.N)
  println("num_subpops: ", sr.num_subpops)
  println("ne: ", sr.ne)
  println("num_attributes: ", sr.num_attributes)
  println("normal_stddev: ", sr.normal_stddev)
  println("ngens: ", sr.ngens)
  println("fitness_mean: ", sr.fitness_mean)
  println("fitness_variance: ", sr.fitness_variance)
  println("attiribute_variance: ", sr.attribute_variance)
end

function run_trials( ) 
  println("stream: ",stream)
  println("stddev: ",normal_stddev())
  sr = spatial_result(N,num_subpops_list[1],ne,num_attributes, ngens, normal_stddev() )
  writeheader( stream, num_subpops_list, sr )
  trial = 1
  for num_subpops in num_subpops_list
    sr = spatial_result(N,num_subpops,ne,num_attributes, ngens, normal_stddev() )
    run_trial(sr)
    writerow(stream,trial,sr)
    print_spatial_result( sr )
    trial += 1
  end
end

function run_trial( sr::spatial_result_type )
  mu = 0.0
  copy_err_prob = 1.0
  (sr.fitness_mean, sr.fitness_variance, sr.attribute_variance ) = spatial_simulation( sr.N, sr.num_subpops, mu, copy_err_prob, 
      sr.ngens, burn_in, sr.ne, sr.num_attributes, normal_stddev(), cdfe, idfe )
end

function writeheader( stream::IO, num_subpops_list::Vector{Int64}, sr::spatial_result_type )
  param_strings = [
    "# $(string(Dates.today()))",
    "# N=$(sr.N)",
    "# num_subpops_list=$(num_subpops_list)",
    "# num_attributes=$(sr.num_attributes)",
    "# ngens=$(sr.ngens)",
    "# burn_in=$(burn_in)",
    "# normal_stddev=$(sr.normal_stddev)"]
  write(stream,join(param_strings,"\n"),"\n")
  heads = [
    "num_subpops",
    "subpop_size",
    "num_emigrants",
    "mean_fitness",
    "variance_fitness",
    "attribute_variance"]
  write(stream,join(heads,","),"\n")
end
    
function writerow( stream::IO, trial::Int64, sr::spatial_result_type )
  line = Any[sr.num_subpops,
          Int(ceil(sr.N/sr.num_subpops)),
          sr.ne,
          sr.fitness_mean,
          sr.fitness_variance,
          sr.attribute_variance]
  write(stream,join(line,","),"\n")
end

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


if isdefined(:simtype)
#  (avg_means,avg_variances, avg_attr_variances) = spatial_simulation( N, m, mu, copy_err_prob, ngens, burn_in, ne, num_attributes, cdfe, idfe, vtbl, subpop_properties )
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
