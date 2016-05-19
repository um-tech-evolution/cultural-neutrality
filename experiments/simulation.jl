#module NeutCultEvo
# Suggested usage:  julia -p 8 -L simulation.jl run_simulation.jl configs/example
include("../src/NeutralCulturalEvolution.jl")
@everywhere using NeutralCulturalEvolution

#=  Moved to the file  run_simulation.jl
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
    "# N=$(N)",
    "# mu_list=$(mu_list)",
    "# ngens=$(ngens)",
    "# burn_in=$(burn_in)",
  ], "\n"), "\n")
  line = join([
    "trial",
    "mu",
    "prob",
    "theta",
  ], ",")
  write(stream, line, "\n")
end

function writerow(stream, trial, mu, prob, theta )
  #println("trial: ",trial)
  line = join(Any[
    trial,
    mu,
    prob,
    theta
  ], ",")
  write(stream, line, "\n")
  #println("line: ",line)
end

function run_trial( N, mu, ngens, burn_in )
  r = ewens_montecarlo(Int32(100000),pop_counts(neutral_poplist(N,mu,ngens, burn_in=burn_in )[ngens]))
  (mu, r.probability, r.theta_estimate)
end

function run_simulation(simname::AbstractString, T, N, mu_list, ngens, burn_in )
  #uprogress = PM.Progress(T, 1, "Running...", 40)
  stream = open("$(simname).csv", "w")
  writeheader(stream)
  trial_list = Float64[]
  for mu in mu_list
    trial_list = vcat(trial_list,fill(mu,T))
  end
  trials = pmap(mu->run_trial( N, mu, ngens, burn_in ), trial_list )
  #println("trials: ",trials)
  t = 1
  for trial in trials
    writerow(stream,t, trial[1],trial[2],trial[3] )
    flush(stream)
    t += 1
  end
  close(stream)
end

#run_simulation(simname, T, N, mu_list, ngens, burn_in )

#end #module
