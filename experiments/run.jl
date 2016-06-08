#=
Command-line program to run the simulation.
 At this time, runs only on a Linux system due to the dependence on the shared library slatkin.so.
  Note:  The location of the shared library file  slatkin.so  must be in LD_LIBRARY_PATH
 The source code for compiling this shared library is at https://github.com/um-tech-evolution/slatkin-exact-tools
 Suggested usage:  julia -p 8 -L simulation.jl run.jl configs/example1
 Suggested usage:  julia -p 8 -L simulation.jl run.jl configs/example2
 Suggested usage:  julia -p 8 -L simulation.jl run.jl configs/example3
 See analysis.jl  for analysis of the resulting CSV file.
 Suggested usage:  julia analysis.jl configs/example1
=#

if length(ARGS) == 0
  simname = "../experiments/configs/example"
else
  simname = ARGS[1]
end

include("$(simname).jl")


try 
  run_simulation(simname, simtype, T, N_list, mu_list, ngens, burn_in, cpower_list=cpower_list )
catch
  run_simulation(simname, simtype, T, N_list, mu_list, ngens, burn_in )
end
