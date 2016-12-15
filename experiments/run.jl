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
  simname = "../experiments/configs/example0"
else
  simname = ARGS[1]
end

include("$(simname).jl")

if  simtype == 0   # Neutral---no conformity
  run_simulation(simname, simtype, T, n_list, N_mu_list, ngens, burn_in, 
      popsize_multiplier=popsize_multiplier, slat_reps=slat_reps )
elseif  simtype == 1   # power conformity
  run_simulation(simname, simtype, T, n_list, N_mu_list, ngens, burn_in, 
      cprob_list=cprob_list, acprob_list=acprob_list, 
      cpower_list=cpower_list, acpower_list=acpower_list, 
      popsize_multiplier=popsize_multiplier, slat_reps=slat_reps )
elseif  simtype == 2   # Acerbi conformity
  run_simulation(simname, simtype, T, n_list, N_mu_list, ngens, burn_in, 
      cprob_list=cprob_list, acprob_list=acprob_list, 
      acer_flag_list=acer_flag_list, bottom_list=bottom_list, topsize_list=topsize_list, bottomsize_list=bottomsize_list,
      popsize_multiplier=popsize_multiplier, slat_reps=slat_reps )
elseif  simtype == 3   # Nearly Neutral Power Conformist
  run_simulation(simname, simtype, T, n_list, N_mu_list, ngens, burn_in, 
      cprob_list=cprob_list, acprob_list=acprob_list, 
      cpower_list=cpower_list, acpower_list=acpower_list, 
      dfe=dfe, popsize_multiplier=popsize_multiplier, slat_reps=slat_reps )
else
  error("illegal simtype")
end
