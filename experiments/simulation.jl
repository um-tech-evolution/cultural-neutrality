#= Nearly neutral version of simulation
=#

#module NeutCultEvo # Should this be a module?
#=
At this time, runs only on a Linux system due to the dependence on the shared library slatkin.so.
 The source code for compiling this shared library is at https://github.com/um-tech-evolution/slatkin-exact-tools
 Suggested usage:  julia -p 8 -L simulation.jl run.jl configs/example1
 Suggested usage:  julia -p 8 -L simulation.jl run.jl configs/example2
 Suggested usage:  julia -p 8 -L simulation.jl run.jl configs/example3
 Note:  These write CSV files such as configs/example1.csv.
 See analysis.jl  for interpretation of the resulting CSV file.
 Suggested usage:  julia analysis.jl configs/example1
=#
include("../src/NeutralCulturalEvolution.jl")
include("conformist.jl")
#Note:  The location of the shared library file  slatkin.so  must be in LD_LIBRARY_PATH on a Linux system
@everywhere using NeutralCulturalEvolution

# If the next 2 lines are changed, the corresponding changes also need to be made in the trial_result type and in writerow.
const p_homoz_headers = [ "p_1_4", "p_1_6", "p_2_6", "p_3_0" ]
const p_homoz_coeffs  = [ 1.4, 1.6, 2.6, 3.0 ]

@doc """ function writeheader(stream::IO, simtype::Int64, T::Int64, n_list::Vector{Int64}, N_mu_list::Vector{Float64}, ngens::Int64, 
    burn_in::Float64, slat_reps::Int64=100000 ) 
Write header line for output CSV file that corresponds to stream.
See the comments for run_simulation for a description of the parameters.
"""
function writeheader(stream::IO, simtype::Int64, T::Int64, n_list::Vector{Int64}, N_mu_list::Vector{Float64}, ngens::Int64, 
    burn_in::Float64, slat_reps::Int64=100000 ) 
  param_strings = [
    "# $(string(Dates.today()))",
    "# trials=$(T)",
    "# n_list=\"$(n_list)\"",
    "# N_mu_list=\"$(N_mu_list)\"",
    "# ngens=$(ngens)",
    "# burn_in=$(burn_in)",
    "# slat_reps=$(slat_reps)",
    "# popsize_multiplier=$(popsize_multiplier)"
    ]
  write(stream, join(param_strings, "\n"), "\n")
  first_heads = ["trial", "N", "n", "N_mu","K"]
  last_heads = vcat([ "s_prob", 
    #"se_prob", 
    #"we_prob", 
    "s_homoz", "w_homoz"], p_homoz_headers )
  if simtype == 1
    mid_heads = ["cpower"]
  elseif simtype == 2
    mid_heads = ["acer_C", "acer_topsize"]
  elseif simtype == 3
    mid_heads = ["nn_select"]
  else
    mid_heads = []
  end
  line = join(vcat( first_heads, mid_heads, last_heads), ",")
  write(stream, line, "\n")
end

# trial result type that includes Watterson homoz, Slatkin homoz, Slatkin prob, and Ewens K estimate
type trial_result
  N::Int64
  n::Int64
  N_mu::Float64
  K::Int64
  cpower::Float64
  acer_C::Float64
  acer_topsize::Int64
  nn_selection::Int64
  s_prob::Float64
  # The following 2 lines can be uncommented if 2 correponding lines uncommented above, and if some lines in hyptest.jl are uncommented.
  #se_prob::Float64
  #we_prob::Float64
  w_homoz::Float64
  s_homoz::Float64
  p_1_4      # p-homozygosity for p = 1.4
  p_1_6      # p-homozygosity for p = 1.6
  p_2_6      # p-homozygosity for p = 2.6
  p_3_0      # p-homozygosity for p = 3.0
end

@doc """ function writerow(stream::IOStream, trial::Int64, N::Int64, N_mu::Float64, K::Int64, theta::Float64 )
Write a row to the CSV file.
"""
function writerow(stream::IO, simtype::Int64, trial::Int64, tr::trial_result  )
  first = Any[
    trial,
    tr.N,           # popsize
    tr.n,           # sample size
    tr.N_mu,        # N/mu
    tr.K,           # The number of alleles in the sample
  ]
  last = Any[
    tr.s_prob,      # Slatkin estimated probability
    # The following 2 lines can be uncommented if 2 correponding lines uncommented above, and if some lines in hyptest.jl are uncommented.
    #tr.se_prob,     # Slatkin exact probability
    #tr.we_prob,     # Watterson exact probability
    tr.s_homoz,     # Slatkin estimated homozygosity = 1/(1+theta) where theta = Slatkin est theta
    tr.w_homoz,     # Watterson estimated homozygosity
    tr.p_1_4,      # p-homozygosity for p = 1.4
    tr.p_1_6,      # p-homozygosity for p = 1.6
    tr.p_2_6,      # p-homozygosity for p = 2.6
    tr.p_3_0,      # p-homozygosity for p = 3.0
  ]
  if simtype == 0 
    mid = Any[]
  elseif simtype == 1
    mid = Any[ 
      tr.cpower 
    ]
  elseif simtype == 2
    mid = Any[ 
      tr.acer_C,
      tr.acer_topsize
    ]
  elseif simtype == 3
    mid = Any[ 
      tr.nn_selection 
    ]
  end
  line = join( vcat( first, mid, last ), "," )
  write(stream, line, "\n")
end

@doc """ function run_trial( n::Int64, N_mu::Float64, ngens::Int64; psize_mult::Float64=1.0, burn_in::Float64=2.0, 
      slat_reps::Int64=10000, cpower::Float64=0.0, acer_C::Float64=0.0, acer_topsize::Int64=10 )
Run a population evolution and return various statitics based on a single generation after burn in.
"""
function run_trial( n::Int64, N_mu::Float64, ngens::Int64; psize_mult::Float64=1.0, burn_in::Float64=2.0, 
      slat_reps::Int64=10000, cpower::Float64=0.0, acer_C::Float64=0.0, acer_topsize::Int64=10, nn_select::Int64=0 )
  if psize_mult > 1.0  # Define popsize N different from sample size n
    N = Int(ceil(psize_mult*n))
    #println("N: ", N)
  else
    N = n
  end
  if cpower==0.0 && acer_C==0.0 && !isdefined(:dfe)
    pop = neutral_poplist(N, N_mu/N, ngens, burn_in=burn_in )[ngens]
  elseif cpower != 0.0
    pop = power_conformist_poplist(N, N_mu/N, ngens, burn_in=burn_in, conformist_power=cpower )[ngens]
  elseif acer_C != 0.0
    pop = acerbi_conformist_poplist(N, N_mu/N, ngens, acer_C, burn_in=burn_in, toplist_size=acer_topsize  )[ngens]
  elseif isdefined(:dfe)   # Nearly neutral
    pop = nearly_neutral_poplist(N, N_mu/N, ngens, dfe, burn_in=burn_in, nn_select=nn_select )[ngens]
  end
  if psize_mult > 1.0 # take a random sample of size n from population of size N
    new_pop = sample_population( pop, n )
    #println("sampling population with sample size: ", n)
  else
    new_pop = pop
  end
  p64 = pop_counts64(new_pop)
  p32 = pop_counts32(new_pop)
  p8 = pop_counts8(new_pop)
  sr = ewens_montecarlo(Int32(slat_reps),p32)
  #se_prob = slatkin_exact( p8, Btbl )
  #we_prob = watterson_exact( p8, Btbl )
  s_homoz = 1.0/(1.0+sr.theta_estimate)
  w_homoz = watterson_homozygosity(p64)
  #ewens_est_K = ewens_K_est( 2*N_mu, N )
  #slatkin_est_K = ewens_K_est( sr.theta_estimate, N )
  #TODO:  include the p_homozygosity parameters in the initialization data
  #p_1_4 = p_homozygosity( p64, 1.4 )
  #p_1_6 = p_homozygosity( p64, 1.6 )
  #p_2_6 = p_homozygosity( p64, 2.6 )
  #p_3_0 = p_homozygosity( p64, 3.0 )
  #trial_result(N, n, N_mu, length(p64), cpower, acer_C, acer_topsize, sr.probability, w_homoz, s_homoz, p_1_4, p_1_6, p_2_1, p_2_4 )
  trial_result(N, n, N_mu, length(p64), cpower, acer_C, acer_topsize, nn_select, sr.probability, 
    # The following 2 lines can be uncommented if 2 correponding lines uncommented above, and if some lines in hyptest.jl are uncommented.
    #se_prob, 
    #we_prob, 
    w_homoz, s_homoz, 
    p_homozygosity( p64, p_homoz_coeffs[1]), 
    p_homozygosity( p64, p_homoz_coeffs[2]), 
    p_homozygosity( p64, p_homoz_coeffs[3]), 
    p_homozygosity( p64, p_homoz_coeffs[4]) )

end

type neutral_trial_type
  n::Int64
  N_mu::Float64
end

type nearly_neutral_trial_type
  n::Int64
  N_mu::Float64
  nn_select::Int64   # 1 means there is nearly neutral selection, 0 means no selection
end

type pconform_trial_type
  n::Int64
  N_mu::Float64
  cpower::Float64
end

type aconform_trial_type
  n::Int64
  N_mu::Float64
  acer_C::Float64
  acer_topsize::Int64
end

function build_global_vars( n_list::Vector{Int64} )
  @everywhere max_n = maximum(n_list)
  @everywhere Btbl = BT(max_n,max_n)
end
  
@doc """ function run_simulation(simname::AbstractString, simtype::Int64, T::Int64, n_list::Vector{Int64}, N_mu_list::Vector{Float64}, 
    ngens::Int64, burn_in::Float64; cpower_list::Vector{Float64}=Float64[], slat_reps::Int64=100000, acer_C_list::Vector{Float64}=Float64[], 

This is the function that actually runs the simulation.  
A list of "trials" (arguments to run_trial()) is produced that is fed to pmap so that trials are run in parallel.
Arguments:
simtype:  1 means use power model of conformity (corresponds to length(cpower_list)>0)
          2 means use Acerbi model of conformity (corresponds to length(acer_C_list)>0)
T:        the number of trials to run for each setting of N and N_mu
n_list:   a list of sample sizes
N_mu_list:  a list of per-generation mutation rates (the per-locus mutation rate is N_mu/N)
ngens:    the number ogenerations to run (not including burn-in)
burn_in:  the number of preliminary generations run to stabilize as a multiple of pop size N
ngens:    the number of generations to run after burn in
slat_reps: the number of monte-carlo reps to use in the Slatkin test
"""

function run_simulation(simname::AbstractString, simtype::Int64, T::Int64, n_list::Vector{Int64}, N_mu_list::Vector{Float64}, 
    ngens::Int64, burn_in::Float64; cpower_list::Vector{Float64}=Float64[], slat_reps::Int64=100000, acer_C_list::Vector{Float64}=Float64[], 
    acer_topsize::Int64=10, popsize_multiplier::Float64=1.0 )
  for n in n_list
    if n > typemax(ConfigInt)
      error("Reset the type ConfigInt in src/aliases.jl so that n <= typemax(ConfigInt) for all n in n_list.")
    end
  end
  build_global_vars( n_list )
  stream = open("$(simname).csv", "w")
  if popsize_multiplier > 1.0
    N_list = map(x->popsize_multiplier*x, n_list)
    println(" N_list: ",N_list)
  end
  writeheader(stream,simtype,T,n_list,N_mu_list,ngens,burn_in)
  if simtype == 0   # neutral, no conformity
    t_list = neutral_trial_type[]
    for n in n_list
      for N_mu in N_mu_list
          push!(t_list,neutral_trial_type(n,N_mu))
      end
    end
  elseif simtype == 1   # power conformity
    t_list = pconform_trial_type[]
    for n in n_list
      for N_mu in N_mu_list
        for cpower in cpower_list
          push!(t_list,pconform_trial_type(n,N_mu,cpower))
        end
      end
    end
  elseif simtype == 2   # Acerbi conformity
    t_list = aconform_trial_type[]
    for n in n_list
      for N_mu in N_mu_list
        for acer_C in acer_C_list
          push!(t_list,aconform_trial_type(n,N_mu,acer_C,acer_topsize))
        end
      end
    end
  elseif simtype == 3   # Nearly neutral
    t_list = nearly_neutral_trial_type[]
    for n in n_list
      for N_mu in N_mu_list
        for nn_sel = 0:1   # 0 means no selection, 1 means nearly neutral selection
          push!(t_list,nearly_neutral_trial_type(n,N_mu,nn_sel))
        end
      end
    end
  end
  trial_list = t_list
  for t in 1:(T-1)
    trial_list = vcat(trial_list,t_list)
  end
  if simtype == 0   # neutral conform
    trial_results = pmap(tr->run_trial( tr.n, tr.N_mu, ngens, burn_in=burn_in, psize_mult=popsize_multiplier, slat_reps=slat_reps ), trial_list )
  elseif simtype == 1   # power conform
    trial_results = pmap(tr->run_trial( tr.n, tr.N_mu, ngens, burn_in=burn_in, psize_mult=popsize_multiplier, slat_reps=slat_reps, cpower=tr.cpower ), trial_list )
  elseif simtype == 2  # Acerbi conform
    trial_results = pmap(tr->run_trial( tr.n, tr.N_mu, ngens, burn_in=burn_in, psize_mult=popsize_multiplier, slat_reps=slat_reps, acer_C=tr.acer_C, acer_topsize=tr.acer_topsize ), trial_list )
  elseif simtype == 3  # Nearly neutral
    trial_results = pmap(tr->run_trial( tr.n, tr.N_mu, ngens, burn_in=burn_in, psize_mult=popsize_multiplier, slat_reps=slat_reps, nn_select=tr.nn_select ), trial_list )
  end
  trial = 1
  for t_result in trial_results
    writerow(stream, simtype, trial, t_result )
    flush(stream)
    trial += 1
  end
  close(stream)
end

#end #module
