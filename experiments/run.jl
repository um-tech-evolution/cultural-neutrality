#=
Command-line program to run the simulation.
# Suggested usage:  julia -p 8 -L simulation.jl run_simulation.jl configs/example
Note:  The location of the shared library file  slatkin.so  must be in LD_LIBRARY_PATH
=#
if length(ARGS) == 0
  simname = "../experiments/configs/example"
else
  simname = ARGS[1]
end

include("$(simname).jl")

run_simulation(simname, simtype, N_list, mu_list, ngens, burn_in )
