include("../src/aliases.jl")
include("../src/freq_scaled_fitness.jl")
include("../src/neutral_poplist.jl")
include("../src/innovation.jl")
include("../src/innovation_collection.jl")
include("../src/nearly_neutral_poplist.jl")
if length(ARGS) == 0
  simname = "../experiments/nn_configs/nn_example1"
else
  simname = ARGS[1]
end
include("$(simname).jl")
println("simname: ",simname)
stream = open("$(simname).csv","w")
println("stream: ",stream)

current_dir = pwd()
date_string = "../data/"*Dates.format(now(),"mm_dd_yy")*"/"
try mkdir(date_string) catch end   # create today's directory with no error if it already exists
println("date string: ",date_string)


@doc """ type trial_result
  A universal trial result for all run_trials functions
  Constructors for specific cases are defined below.
"""
type trial_result
  nn_simtype::Int64
  N::Int64
  N_mu::Float64
  L::Int64    # number of loci
  ngens::Int64
  burn_in::Float64
  fix_minimum::Float64
  dfe::Function
  dfe_str::AbstractString
  number_active::Int64
  number_extinctions::Int64
  number_fixations::Int64
  extinction_time::Float64
  fixation_time::Float64
  average_fitness_fixed::Float64
  average_fitness_extinct::Float64
  average_fitness_all::Float64
end

# Constructor that sets the parameters
function trial_result( nn_symtype::Int64, N::Int64, N_mu::Float64, L::Int64, ngens::Int64, 
    fix_minimum::Float64, burn_in::Float64, dfe::Function, dfe_str::AbstractString )
  trial_result( nn_simtype, N, N_mu, L, ngens, burn_in, fix_minimum, dfe, dfe_str, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0 )
end

function print_trial_result( tr::trial_result )
  if tr.nn_simtype == 1
    println("\ninfinite alleleles model")
  else
    println("nn_simtype: ", tr.nn_simtype)
  end
  println("N: ", tr.N)
  println("N_mu: ", tr.N_mu)
  println("L: ", tr.L)
  println("ngens: ", tr.ngens)
  println("burn_in: ", tr.burn_in)
  println("fix_minimum: ", tr.fix_minimum)
  println("dfe: ", tr.dfe)
  println("dfe_str: ", tr.dfe_str)
  println("number_active: ", tr.number_active)
  println("number_extinctions: ", tr.number_extinctions)
  println("number_fixations: ", tr.number_fixations)
  println("extinction_time: ", tr.extinction_time)
  println("fixation_time: ", tr.fixation_time)
  println("average_fitness_fixed: ", tr.average_fitness_fixed)
  println("average_fitness_extinct: ", tr.average_fitness_extinct)
  println("average_fitness_all: ", tr.average_fitness_all)
end

# Convert an innovation_oollection to a trial_result
# The innovation_connection and the trial_result must agree in fields n, N, 
function convert_to_trial_result( ic::innovation_collection, tr::trial_result )
  tr.number_active = length(ic.active)
  tr.number_extinctions = length(ic.extinct)
  tr.number_fixations = length(ic.fixed)
  tr.fix_minimum = ic.fix_minimum
  tr.extinction_time = average_time_to_extinction( ic )
  tr.fixation_time = average_time_to_fixation( ic )
  tr.average_fitness_fixed = average_fitness_fixed( ic::innovation_collection )
  tr.average_fitness_extinct = average_fitness_extinct( ic::innovation_collection )
  tr.average_fitness_all = average_fitness_all( ic::innovation_collection )
end

function run_trials()
  println("stream: ",stream)
  trial = 1
  writeheader(stream, N_list, N_mu_list, ngens, burn_in )
  for N_mu in N_mu_list
    for N in N_list
      tr = trial_result( nn_simtype, N, N_mu, L, ngens, fix_minimum, burn_in, dfe, dfe_str )
      run_trial( tr )
      writerow(stream, trial, tr )
      trial += 1
    end
  end
end

function run_trial( tr::trial_result )
  if tr.nn_simtype == 1
    ic = innovation_collection( tr.fix_minimum )
    poplist = nearly_neutral_poplist(tr.N,tr.N_mu,tr.ngens,tr.dfe,combine=false,ic=ic)
    convert_to_trial_result( ic, tr )
    print_trial_result( tr )
    return tr
  else
    println("nn_simtype ",nn_simtype," not implemented")
  end
end

@doc """ function writeheader()
    burn_in::Float64, slat_reps::Int64=100000 ) 
Write header line for output CSV file that corresponds to stream.
See the comments for run_simulation for a description of the parameters.
"""
function writeheader(stream::IO, N_list::Vector{Int64}, N_mu_list::Vector{Float64}, ngens::Int64, 
    burn_in::Float64) 
  dfe_params = [:dfe_adv_prob, :dfe_adv_alpha, :dfe_adv_theta, :dfe_disadv_prob, :dfe_disadv_alpha, :dfe_disadv_theta]
  param_strings = [
    "# $(string(Dates.today()))",
    "# $((nn_simtype==1)?"infinite alleles model":"infinite_sites_model")",
    "# N_list=\"$(N_list)\"",
    "# N_mu_list=\"$(N_mu_list)\"",
    "# ngens=$(ngens)",
    "# burn_in=$(burn_in)",
    "# fix_minimum=$(fix_minimum)",
    "# dfe=$(dfe)",
    "# dfe_str=$(dfe_str)"
    ]
  for dp in dfe_params
    if isdefined(dp)
      Base.push!(param_strings,"# $(string(dp))=$(eval(dp))")
    end
  end
  write(stream, join(param_strings, "\n"), "\n")
  first_heads = ["trial", "N", "N_mu","ngens"]
  mid_heads = []
  last_heads =
  [ "num_extinct", 
    "num_fixed", 
    "fraction_fixed",
    "ave_extinct_time", 
    "ave_fixed_time", 
    "ave_extinct_selcoef", 
    "ave_fixed_selcoef", 
    "ave_all_selcoef"] 
  #=
  if tr.nn_simtype == 0  # Neutral
    mid_heads = []
  elseif tr.nn_simtype == 1  # Power conformist
    mid_heads = [ "cprob", "acprob", "cpower", "acpower"]
  else   
    mid_heads = []
  end
  =#
  line = join(vcat( first_heads, mid_heads, last_heads), ",")
  write(stream, line, "\n")
end

function writerow(stream::IO, trial::Int64, tr::trial_result  )
  first = Any[
    trial,
    tr.N,           # popsize
    tr.N_mu,        # N/mu, population mutation rate
    tr.ngens,       # The number of generations after burn in
  ]
  if tr.nn_simtype == 0
    mid = Any[]
  elseif tr.nn_simtype == 1
    mid = Any[
      tr.number_extinctions,
      tr.number_fixations,
      Float64(tr.number_fixations)/(tr.number_extinctions+tr.number_fixations),
      tr.extinction_time,
      tr.fixation_time,
      tr.average_fitness_extinct,
      tr.average_fitness_fixed,
      tr.average_fitness_all
    ]
  end
  line = join( vcat( first, mid ), "," )
  write(stream, line, "\n")
end


try
  run_trials()
  close(stream)
except
  close(stream)
end
