#module NeutCultEvo
# Suggested usage:  julia -p 8 -L simulation.jl run_simulation.jl configs/example
include("../src/NeutralCulturalEvolution.jl")
@everywhere using NeutralCulturalEvolution

#=  Moved to the file:  run_simulation.jl
if length(ARGS) == 0
  simname = "../experiments/configs/example"
else
  simname = ARGS[1]
end

@everywhere include("$(simname).jl")
=#

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
    "theta",
  ], ",")
  write(stream, line, "\n")
end

function writerow(stream, trial, N, mu, prob, theta )
  #println("trial: ",trial)
  line = join(Any[
    trial,
    N,
    mu,
    prob,
    theta
  ], ",")
  write(stream, line, "\n")
  #println("line: ",line)
end

function run_trial_slatkin( N, mu, ngens, burn_in, slat_reps )
  r = ewens_montecarlo(Int32(slat_reps),pop_counts32(neutral_poplist(N,mu/N,ngens, burn_in=burn_in )[ngens]))
  (N, mu, r.probability, r.theta_estimate)
end

function run_trial_watterson( N, mu, ngens, burn_in )
  theta_estimate = watterson(pop_counts64(neutral_poplist(N,mu/N,ngens, burn_in=burn_in )[ngens]))
  (N, mu, 0, theta_estimate)
end

function run_simulation(simname::AbstractString, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, 
    ngens::Int64, burn_in::Float64, slat_reps::Int64 )
  #uprogress = PM.Progress(T, 1, "Running...", 40)
  #println("rs burn_in: ",burn_in)
  stream = open("$(simname).csv", "w")
  writeheader(stream)
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
  #println("trial list: ",trial_list)
  if slat_reps == 0
    trials = pmap(tr->run_trial_watterson( tr[1], tr[2], ngens, burn_in ), trial_list )
  else
    trials = pmap(tr->run_trial_slatkin( tr[1], tr[2], ngens, burn_in, slat_reps ), trial_list )
  end
  #println("trials: ",trials)
  t = 1
  for trial in trials
    writerow(stream, t, trial[1],trial[2],trial[3],trial[4] ) #trial[1]:N,  trial[2]:mu, trial[3]:r.probability
    flush(stream)
    t += 1
  end
  close(stream)
end

#run_simulation(simname, T, N, mu_list, ngens, burn_in )

#end #module
