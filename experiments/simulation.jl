#module NeutCultEvo
include("../src/NeutralCulturalEvolution.jl")
import ProgressMeter
const PM = ProgressMeter
import NeutralCulturalEvolution

if length(ARGS) == 0
  simname = "configs/example"
else
  simname = ARGS[1]
end

include("$(simname).jl")

function writeheader(stream)
  write(stream, join([
    "# trials=$(T)",
    "# N=$(N)",
    "# mu=$(mu)",
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

function runtrial(trial, stream, progress )
  plist = neutral_poplist( N, mu, ngens, burn_in=burn_in )
  result = ewens_montecarlo(100000, pop_counts(plist[ngens]) )
  writerow(stream, trial, result.probability, result.theta_estimate )
  PM.next!(progress)
end

function run_simulation(simname::AbstractString)
  progress = PM.Progress(T, 1, "Running...", 40)
  stream = open("$(simname).csv", "w")
  writeheader(stream)
  for trial = 1:T
    runtrial(trial, stream, progress )
    flush(stream)
  end
  close(stream)
end

run_simulation(simname)

#end #module
