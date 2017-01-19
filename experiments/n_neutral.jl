include("../src/NeutralCulturalEvolution.jl")
#=
include("../src/aliases.jl")
include("../src/freq_scaled_fitness.jl")
include("../src/neutral_poplist.jl")
include("../src/innovation.jl")
include("../src/innovation_collection.jl")
include("../src/nearly_neutral_poplist.jl")
=#
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
  expected_richness::Float64  # sum_{i=0}^{n-1} theta/(theta+i) where theta = 2*N_mu.
  average_richness::Float64   # Average number of traits in populations
  w_homoz::Float64
  IQV::Float64
  number_active::Int64
  number_extinctions::Int64
  number_fixations::Int64
  extinction_time::Float64
  fixation_time::Float64
  average_fitness_fixed::Float64
  average_fitness_extinct::Float64
  average_fitness_all::Float64
  t2::Float64  # Average turnover for y=2
  t4::Float64  # Average turnover for y=4
  t8::Float64  # Average turnover for y=8
  t12::Float64  # Average turnover for y=12
  t20::Float64  # Average turnover for y=20
end

function add_expected_richness( tr::trial_result )
  sum = 0.0
  theta = 2*tr.N_mu
  for i = 0:(tr.N-1)
    sum += 1.0/(theta+i)
  end
  tr.expected_richness = theta*sum
end

# Constructor that sets the parameters
function trial_result( nn_symtype::Int64, N::Int64, N_mu::Float64, L::Int64, ngens::Int64,
    fix_minimum::Float64=0.5, burn_in::Float64=2.0, dfe::Function=dfe_neutral, 
    dfe_str::AbstractString="neutral" )
  tr = trial_result( nn_simtype, N, N_mu, L, ngens, burn_in, fix_minimum, dfe, dfe_str, 0.0,
      0.0, 0.0, 0.0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0,0,0.0 )
  add_expected_richness( tr )
  tr
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
  println("expected richness: ",tr.expected_richness)
  println("average richness: ",tr.average_richness)
  println("w_homoz: ",tr.w_homoz)
  println("IQV: ",tr.IQV)
  println("number_active: ", tr.number_active)
  println("number_extinctions: ", tr.number_extinctions)
  println("number_fixations: ", tr.number_fixations)
  println("extinction_time: ", tr.extinction_time)
  println("fixation_time: ", tr.fixation_time)
  println("average_fitness_fixed: ", tr.average_fitness_fixed)
  println("average_fitness_extinct: ", tr.average_fitness_extinct)
  println("average_fitness_all: ", tr.average_fitness_all)
  println("turnover 2: ",tr.t2)
  println("turnover 4: ",tr.t4)
  println("turnover 8: ",tr.t8)
  println("turnover 12: ",tr.t12)
  println("turnover 20: ",tr.t20)
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
  tr = trial_result( nn_simtype, N_list[1], N_mu_list[1], L, ngens, fix_minimum, burn_in, dfe, dfe_str )
  writeheader(stream, N_list, N_mu_list, tr )
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
    add_stats_to_trial_result!( tr, poplist )
    print_trial_result( tr )
    return tr
  else
    println("nn_simtype ",nn_simtype," not implemented")
  end
end

function add_stats_to_trial_result!( tr::trial_result, poplist::Vector{Population} )
  pcounts = map(pop_counts64,poplist)
  tr.average_richness = mean(map(x->length(x),pcounts))
  tr.w_homoz = mean(map(watterson_homozygosity,pcounts))
  tr.IQV = mean(map(IQV,pcounts))
  turnover_list = turnover(poplist,[3,6,10,20,40],tr.N_mu)
  len = minimum(map(length,turnover_list))
  sum_list = zeros(Int64,len)
  for t in turnover_list
    for i = 1:len
      sum_list[i] += t[i]
    end
  end
  tr.t2 = len >= 1 ? Float64(sum_list[1])/length(turnover_list) : -1.0
  tr.t4 = len >= 2 ? Float64(sum_list[2])/length(turnover_list) : -1.0
  tr.t8 = len >= 2 ? Float64(sum_list[2])/length(turnover_list) : -1.0
  tr.t12 = len >= 3 ? Float64(sum_list[3])/length(turnover_list) : -1.0
  tr.t20 = len >= 4 ? Float64(sum_list[4])/length(turnover_list) : -1.0
end

@doc """ function writeheader()
    burn_in::Float64, slat_reps::Int64=100000 ) 
Write header line for output CSV file that corresponds to stream.
See the comments for run_simulation for a description of the parameters.
"""
function writeheader(stream::IO, N_list::Vector{Int64}, N_mu_list::Vector{Float64}, tr::trial_result )
  dfe_params = [:dfe_adv_prob, :dfe_adv_alpha, :dfe_adv_theta, :dfe_disadv_prob, :dfe_disadv_alpha, :dfe_disadv_theta]
  param_strings = [
    "# $(string(Dates.today()))",
    "# $((nn_simtype==1)?"infinite alleles model":"infinite_sites_model")",
    "# N_list=\"$(N_list)\"",
    "# N_mu_list=\"$(N_mu_list)\"",
    "# ngens=$(tr.ngens)",
    "# burn_in=$(tr.burn_in)",
    "# fix_minimum=$(tr.fix_minimum)",
    "# dfe=$(tr.dfe)",
    "# dfe_str=$(tr.dfe_str)"
    ]
  for dp in dfe_params
    if isdefined(dp)
      Base.push!(param_strings,"# $(string(dp))=$(eval(dp))")
    end
  end
  write(stream, join(param_strings, "\n"), "\n")
  first_heads = ["trial", "N", "N_mu"]
  mid_heads = []
  last_heads =
  [ "expected_richness",
    "average_richness",
    "w_homoz",
    "IQV",
    "num_extinct", 
    "num_fixed", 
    "fraction_fixed",
    "ave_extinct_time", 
    "ave_fixed_time", 
    "ave_extinct_selcoef", 
    "ave_fixed_selcoef", 
    "ave_all_selcoef",
    "t2",
    "t4",
    "t8",
    "t12",
    "t20"
  ] 
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
  ]
  if tr.nn_simtype == 0
    mid = Any[]
  elseif tr.nn_simtype == 1
    mid = Any[
      tr.expected_richness,
      tr.average_richness,
      tr.w_homoz,
      tr.IQV,
      tr.number_extinctions,
      tr.number_fixations,
      Float64(tr.number_fixations)/(tr.number_extinctions+tr.number_fixations),
      tr.extinction_time,
      tr.fixation_time,
      tr.average_fitness_extinct,
      tr.average_fitness_fixed,
      tr.average_fitness_all,
      tr.t2,
      tr.t4,
      tr.t8,
      tr.t12,
      tr.t20
    ]
  end
  line = join( vcat( first, mid ), "," )
  write(stream, line, "\n")
end

run_trials()
#=
try
  run_trials()
  close(stream)
except
  close(stream)
end
=#


