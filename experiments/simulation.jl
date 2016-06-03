# Should this be a module?
#module NeutCultEvo
# Suggested usage:  julia -p 8 -L simulation.jl run_simulation.jl configs/example
include("../src/NeutralCulturalEvolution.jl")
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
function writeheader(stream)
  write(stream, join([
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

# For estimating K
function writeheaderK(stream)
  write(stream, join([
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
    "K",
    "est_theta",
    "true_theta",
  ], ",")
  write(stream, line, "\n")
end

# TODO:  Write documentation
function writerow(stream, trial, N, mu, prob, theta )
  #println("trial: ",trial)
  line = join(Any[
    trial,
    N,
    mu,
    prob,
    theta,     # estimated theta
    2.0*mu     # true theta
  ], ",")
  write(stream, line, "\n")
  #println("line: ",line)
end

# For estimating K
function writerowK(stream, trial, N, mu, K, theta )
  #println("trial: ",trial)
  line = join(Any[
    trial,
    N,
    mu,
    K,
    theta,
    2.0*mu
  ], ",")
  write(stream, line, "\n")
  #println("line: ",line)
end

# TODO:  Write documentation
# TODO:  Change so that a type with names is returned instead of a tuple
function run_trial_slatkin( N, mu, ngens, burn_in, slat_reps )
  r = ewens_montecarlo(Int32(slat_reps),pop_counts32(neutral_poplist(N,mu/N,ngens, burn_in=burn_in )[ngens]))
  (N, mu, r.probability, r.theta_estimate)
end

# TODO:  Write documentation
# TODO:  Change so that a type with names is returned instead of a tuple
function run_trial_watterson( N, mu, ngens, burn_in )
  theta_estimate = watterson(pop_counts64(neutral_poplist(N,mu/N,ngens, burn_in=burn_in )[ngens]))
  (N, mu, 0, theta_estimate)
end

# TODO:  Write documentation
# TODO:  Change so that a type with names is returned instead of a tuple
function run_trial_K( N, mu, ngens, burn_in )
  pop_counts = pop_counts64(neutral_poplist(N,mu/N,ngens, burn_in=burn_in )[ngens])
  theta_estimate = watterson(pop_counts)
  K_estimate = length(pop_counts)
  (N, mu, K_estimate, theta_estimate)
end

@doc """ function run_simulation(simname::AbstractString, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, 
    ngens::Int64, burn_in::Float64, slat_reps::Int64 )

This is the function that actually runs the simulation.  
A list of "trials" (arguments to run_trial()) is produced that is fed to pmap so that trials are run in parallel.

Note that "slat_reps" is used as a switch for the kind of trial:
slat_reps > 0  means Slatkin trials with this number of reps
slat_reps == 0  means Watterson trials
slat_reps < 0  means K estimation trials (K is the number of types of alleles)
"""

function run_simulation(simname::AbstractString, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, 
    ngens::Int64, burn_in::Float64, slat_reps::Int64 )
  #uprogress = PM.Progress(T, 1, "Running...", 40)  # Don't know how to do this with pmap
  stream = open("$(simname).csv", "w")
  if slat_reps < 0  # K_estimate
    writeheaderK(stream)
  else
    writeheader(stream)
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
  if slat_reps == 0  # Watterson test
    trials = pmap(tr->run_trial_watterson( tr[1], tr[2], ngens, burn_in ), trial_list )
  elseif slat_reps < 0  # K_estimate
    trials = pmap(tr->run_trial_K( tr[1], tr[2], ngens, burn_in ), trial_list )
  else   # Slatkin test
    trials = pmap(tr->run_trial_slatkin( tr[1], tr[2], ngens, burn_in, slat_reps ), trial_list )
  end
  t = 1
  for trial in trials
    if slat_reps < 0  # K_estimate
      #  trial[1]:N,  trial[2]:mu, trial[3]:r.probability or K_estimate,  trial[4]: est_theta
      writerowK(stream, t, trial[1],trial[2],trial[3],trial[4] ) 
    else
      writerow(stream, t, trial[1],trial[2],trial[3],trial[4] ) 
    end
    flush(stream)
    t += 1
  end
  close(stream)
end


#end #module
