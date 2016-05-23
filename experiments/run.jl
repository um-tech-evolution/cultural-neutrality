if length(ARGS) == 0
  simname = "../experiments/configs/example"
else
  simname = ARGS[1]
end

include("$(simname).jl")

run_simulation(simname, T, N_list, mu_list, ngens, burn_in, slat_reps )
