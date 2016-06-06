# Should this be a module?
#module NeutCultEvo
# Suggested usage:  julia -p 8 -L simulation.jl run_simulation.jl configs/example
include("../src/NeutralCulturalEvolution.jl")
#Note:  The location of the shared library file  slatkin.so  must be in LD_LIBRARY_PATH
@everywhere using NeutralCulturalEvolution

#=  Moved to the file:  run.jl
if length(ARGS) == 0
  simname = "../experiments/configs/example"
else
  simname = ARGS[1]
end

@everywhere include("$(simname).jl")
=#

# TODO:  Write documentation
@doc """ function writeheader_Slatkin(stream::IOStream, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, ngens::Int64, 
    burn_in::Float64, slat_reps::Int64=100000 )
Write header line for output CSV file that corresponds to stream.
See the comments for run_simulation for a description of the parameters.
"""
function writeheader_Slatkin(stream::IOStream, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, ngens::Int64, 
    burn_in::Float64, slat_reps::Int64=100000 )
  write(stream, join([
    "# Slatkin",
    "# trials=$(T)",
    "# N_list=$(N_list)",
    "# mu_list=$(mu_list)",
    "# ngens=$(ngens)",
    "# burn_in=$(burn_in)",
    "# slat_reps=$(slat_reps)",
  ], "\n"), "\n")
  line = join([
    "trial",
    "N",
    "N_mu",
    "prob",
    "est_theta",
    "true_theta",
  ], ",")
  write(stream, line, "\n")
end

@doc """ function writeheader_Watterson(stream::IOStream, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, ngens::Int64, burn_in::Float64 )
Write header line for output CSV file that corresponds to stream.
See the comments for run_simulation for a description of the parameters.
"""
function writeheader_Watterson(stream::IOStream, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, ngens::Int64, burn_in::Float64 )
  write(stream, join([
    "# Watterson",
    "# trials=$(T)",
    "# N_list=$(N_list)",
    "# mu_list=$(mu_list)",
    "# ngens=$(ngens)",
    "# burn_in=$(burn_in)",
  ], "\n"), "\n")
  line = join([
    "trial",
    "N",
    "N_mu",
    "prob",
    "est_theta",
    "true_theta",
  ], ",")
  write(stream, line, "\n")
end

# For estimating K
@doc """ function writeheader_estK(stream::IOStream, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, ngens::Int64, burn_in::Float64 )
Write header line for output CSV file that corresponds to stream.
See the comments for run_simulation for a description of the parameters.
"""

function writeheader_estK(stream::IOStream, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, ngens::Int64, burn_in::Float64 )
  write(stream, join([
    "# Estimate K",
    "# trials=$(T)",
    "# N_list=$(N_list)",
    "# mu_list=$(mu_list)",
    "# ngens=$(ngens)",
    "# burn_in=$(burn_in)",
  ], "\n"), "\n")
  line = join([
    "trial",
    "N",
    "N_mu",
    "K",
    "est_theta",
    "true_theta",
  ], ",")
  write(stream, line, "\n")
end



@doc """ function writerow(stream::IOStream, trial::Int64, N::Int64, mu::Float64, prob::Float64, theta::Float64 )
Write a row to the CSV file for Slatkin and Watterson.
"""
function writerow(stream::IOStream, trial::Int64, N::Int64, mu::Float64, prob::Float64, theta::Float64 )
  line = join(Any[
    trial,
    N,
    mu,
    prob,
    theta,     # estimated theta
    2.0*mu     # true theta
  ], ",")
  write(stream, line, "\n")
end

# For estimating K
@doc """ function writerow_estK(stream::IOStream, trial::Int64, N::Int64, mu::Float64, K::Int64, theta::Float64 )
Write a row to the CSV file for K estimation.
"""
function writerow_estK(stream::IOStream, trial::Int64, N::Int64, mu::Float64, K::Int64, theta::Float64 )
  line = join(Any[
    trial,
    N,
    mu,
    K,
    theta,
    2.0*mu
  ], ",")
  write(stream, line, "\n")
end

# TODO:  Write documentation
# TODO:  Change so that a type with names is returned instead of a tuple
function run_trial_slatkin( N::Int64, mu::Float64, ngens::Int64, burn_in::Float64, slat_reps::Int64 )
  r = ewens_montecarlo(Int32(slat_reps),pop_counts32(neutral_poplist(N,mu/N,ngens, burn_in=burn_in )[ngens]))
  (N, mu, r.probability, r.theta_estimate)
end

# TODO:  Write documentation
# TODO:  Change so that a type with names is returned instead of a tuple
function run_trial_watterson( N, mu, ngens, burn_in )
  theta_estimate = watterson_theta(pop_counts64(neutral_poplist(N,mu/N,ngens, burn_in=burn_in )[ngens]))
  (N, mu, 0.0, theta_estimate)
end

# TODO:  Write documentation
# TODO:  Change so that a type with names is returned instead of a tuple
function run_trial_K( N::Int64, mu::Float64, ngens::Int64, burn_in::Float64 )
  pop_counts = pop_counts64(neutral_poplist(N,mu/N,ngens, burn_in=burn_in )[ngens])
  theta_estimate = watterson_theta(pop_counts)
  K_estimate = length(pop_counts)
  (N, mu, K_estimate, theta_estimate)
end

@doc """ function run_simulation(simname::AbstractString, simtype::bool, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, 
    ngens::Int64, burn_in::Float64; slat_reps::Int64==10000 )

This is the function that actually runs the simulation.  
A list of "trials" (arguments to run_trial()) is produced that is fed to pmap so that trials are run in parallel.
Arguments:
simtype:  1 means use Slatkin test
          2 means use Watterson test
          3 means estimate K where K is the number of allele types
T:        the number of trials to run for each setting of N and mu
N_list:   a list of population sizes
mu_list:  a list of per-generation mutation rates (the per-locus mutation rate is mu/N)
ngens:    the number ogenerations to run (not including burn-in)
burn_in:  the number of preliminary generations run to stabilize as a multiple of pop size N
ngens:    the number of generations to run after burn in
slat_reps: the number of monte-carlo reps to use in the Slatkin test
"""

function run_simulation(simname::AbstractString, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, 
    ngens::Int64, burn_in::Float64; slat_reps::Int64=100000 )
  #uprogress = PM.Progress(T, 1, "Running...", 40)  # Don't know how to do this with pmap
  stream = open("$(simname).csv", "w")
  if simtype == 3  # K_estimate
    writeheader_estK(stream, T, N_list, mu_list, ngens, burn_in )
  elseif simtype == 2
    writeheader_Watterson(stream, T, N_list, mu_list, ngens, burn_in )
  elseif simtype == 1
    writeheader_Slatkin(stream, T, N_list, mu_list, ngens, burn_in, slat_reps )
  end
  N_mu_list = Tuple{Int64,Float64}[]
  for N in N_list
    for mu in mu_list
      push!(N_mu_list,(N,mu))
    end
  end
  trial_list = N_mu_list
  for t in 1:(T-1)
    trial_list = vcat(trial_list,N_mu_list)
  end
  if simtype == 2  # Watterson test
    trials = pmap(tr->run_trial_watterson( tr[1], tr[2], ngens, burn_in ), trial_list )
  elseif simtype == 3  # K_estimate
    trials = pmap(tr->run_trial_K( tr[1], tr[2], ngens, burn_in ), trial_list )
  elseif simtype == 1   # Slatkin test
    trials = pmap(tr->run_trial_slatkin( tr[1], tr[2], ngens, burn_in, slat_reps ), trial_list )
  end
  t = 1
  for trial in trials
    if simtype == 3  # K_estimate
      #  trial[1]:N,  trial[2]:mu, trial[3]:r.probability or K_estimate,  trial[4]: est_theta
      writerow_estK(stream, t, trial[1],trial[2],trial[3],trial[4] ) 
    else
      writerow(stream, t, trial[1],trial[2],trial[3],trial[4] ) 
    end
    flush(stream)
    t += 1
  end
  close(stream)
end


#end #module
