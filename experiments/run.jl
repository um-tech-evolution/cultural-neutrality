#=
Command-line program to run the simulation.
 At this time, runs only on a Linux system due to the dependence on the shared library slatkin.so.
  Note:  The location of the shared library file  slatkin.so  must be in LD_LIBRARY_PATH
 The source code for compiling this shared library is at https://github.com/um-tech-evolution/slatkin-exact-tools
 Suggested usage:  julia -p 8 -L simulation.jl run.jl configs/example1    # positive power conformity
 Suggested usage:  julia -p 8 -L simulation.jl run.jl configs/example2    # negative power conformity
 Suggested usage:  julia -p 8 -L simulation.jl run.jl configs/example3    # positive Acerbi conformity
 Suggested usage:  julia -p 8 -L simulation.jl run.jl configs/example4    # negative Acerbi conformity
 See hyptestjl  for analysis of the resulting CSV file.
 Suggested usage:  julia hyptest.jl configs/example1
=#

if length(ARGS) == 0
  simname = "../experiments/configs/example1"
else
  simname = ARGS[1]
end

include("$(simname).jl")


if isdefined(:cp_list)
  simtype = 1   # power conformity
  run_simulation(simname, simtype, T, n_list, N_mu_list, ngens, burn_in, cpower_list=cp_list, 
    popsize_multiplier=popsize_multiplier )
elseif isdefined(:acer_C_list) && isdefined(:acer_topsize)
  simtype = 2   # Acerbi conformity
  run_simulation(simname, simtype, T, n_list, N_mu_list, ngens, burn_in, acer_C_list=acer_C_list, 
    acer_topsize=acer_topsize, popsize_multiplier=popsize_multiplier )
elseif isdefined(:dfe) 
  simtype = 3   # Nearly Neutral
  run_simulation(simname, simtype, T, n_list, N_mu_list, ngens, burn_in, popsize_multiplier=popsize_multiplier )
else
  simtype = 0   # Neutral---no conformity
  run_simulation(simname, simtype, T, n_list, N_mu_list, ngens, burn_in, popsize_multiplier=popsize_multiplier )
end
