# Suggested usage:
# julia -p 8 -L simulation.jl run_simulation.jl configs/example
if length(ARGS) == 0
  simname = "../experiments/configs/example"
else
  simname = ARGS[1]
end

include("$(simname).jl")

run_simulation(simname, T, N, mu_list, ngens, burn_in )
