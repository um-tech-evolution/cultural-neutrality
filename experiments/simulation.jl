#module NeutCultEvo # Should this be a module?
#=
At this time, runs only on a Linux system due to the dependence on the shared library slatkin.so.
 The source code for compiling this shared library is at https://github.com/um-tech-evolution/slatkin-exact-tools
 Suggested usage:  julia -p 8 -L simulation.jl run.jl configs/example1
 Suggested usage:  julia -p 8 -L simulation.jl run.jl configs/example2
 Suggested usage:  julia -p 8 -L simulation.jl run.jl configs/example3
 Note:  These write CSV files such as configs/example1.csv.
 See analysis.jl  for interpretation of the resulting CSV file.
 Suggested usage:  julia analysis.jl configs/example1
=#
include("../src/NeutralCulturalEvolution.jl")
#Note:  The location of the shared library file  slatkin.so  must be in LD_LIBRARY_PATH on a Linux system
@everywhere using NeutralCulturalEvolution

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
    "K_est",
    "exp_K",
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
function writerow_estK(stream::IOStream, trial::Int64, N::Int64, mu::Float64, K_est::Float64, theta::Float64 )
  line = join(Any[
    trial,
    N,
    mu,
    K_est,
    ewens_K_est(2.0*mu,N),
    theta,
    2.0*mu
  ], ",")
  write(stream, line, "\n")
end

# trial result type for both Slatkin and Watterson
type sw_trial_result
  N::Int64
  mu::Float64
  prob::Float64
  theta_est::Float64
end

# TODO:  Write documentation
# TODO:  Change so that a type with names is returned instead of a tuple
@doc """ function run_trial_slatkin( N::Int64, mu::Float64, ngens::Int64, burn_in::Float64, slat_reps::Int64 )
Run a neutral population evolution and apply the Slatkin monte carlo test to the result of a single generation.
"""
function run_trial_slatkin( N::Int64, mu::Float64, ngens::Int64, burn_in::Float64, slat_reps::Int64 )
  r = ewens_montecarlo(Int32(slat_reps),pop_counts32(neutral_poplist(N,mu/N,ngens, burn_in=burn_in )[ngens]))
  sw_trial_result(N, mu, r.probability, r.theta_estimate)
end

@doc """ function run_trial_watterson( N, mu, ngens, burn_in )
Run a neutral population evolution and apply the Watterson test to the result of a single generation.
"""
function run_trial_watterson( N, mu, ngens, burn_in )
  theta_estimate = watterson_theta(pop_counts64(neutral_poplist(N,mu/N,ngens, burn_in=burn_in )[ngens]))
  sw_trial_result(N, mu, 0.0, theta_estimate)
end

# trial result type for K estimation
type K_est_trial_result
  N::Int64
  mu::Float64
  K_est::Float64
  theta_est::Float64
end

# TODO:  Write documentation
# TODO:  Change so that a type with names is returned instead of a tuple
@doc """ function run_trial_K( N::Int64, mu::Float64, ngens::Int64, burn_in::Float64 )  # simple_popcounts is a bare-bones version of neutral_poplist
Run a neutral population evolution and apply K estimation to the result of a single generation.
"""
function run_trial_K( N::Int64, mu::Float64, ngens::Int64, burn_in::Float64 )  # simple_popcounts is a bare-bones version of neutral_poplist
  #pop_counts = pop_counts64(simple_poplist(N,mu,ngens, burn_in=burn_in )[ngens])
  pop_counts = pop_counts64(neutral_poplist(N,mu/N,ngens, burn_in=burn_in )[ngens])
  theta_estimate = watterson_theta(pop_counts)
  K_estimate = length(pop_counts)
  K_est_trial_result(N, mu, K_estimate, theta_estimate)
end

type trial_type
  N::Int64
  mu::Float64
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

function run_simulation(simname::AbstractString, simtype::Int64, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, 
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
  #N_mu_list = Tuple{Int64,Float64}[]
  N_mu_list = trial_type[]
  for N in N_list
    for mu in mu_list
      push!(N_mu_list,trial_type(N,mu))
    end
  end
  trial_list = N_mu_list
  for t in 1:(T-1)
    trial_list = vcat(trial_list,N_mu_list)
  end
  if simtype == 2  # Watterson test
    trials = pmap(tr->run_trial_watterson( tr.N, tr.mu, ngens, burn_in ), trial_list )
  elseif simtype == 3  # K_estimate
    trials = pmap(tr->run_trial_K( tr.N, tr.mu, ngens, burn_in ), trial_list )
  elseif simtype == 1   # Slatkin test
    trials = pmap(tr->run_trial_slatkin( tr.N, tr.mu, ngens, burn_in, slat_reps ), trial_list )
  end
  t = 1
  for trial in trials
    if simtype == 3  # K_estimate
      writerow_estK(stream, t, trial.N,trial.mu,trial.K_est,trial.theta_est ) 
    else   # both Slatkin and Watterson tests
      writerow(stream, t, trial.N, trial.mu, trial.prob, trial.theta_est )
    end
    flush(stream)
    t += 1
  end
  close(stream)
end


#end #module
