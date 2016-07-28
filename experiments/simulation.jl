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

@doc """ function writeheader_Slatkin(stream::IOStream, T::Int64, N_list::Vector{Int64}, N_mu_list::Vector{Float64}, ngens::Int64, 
    burn_in::Float64, slat_reps::Int64=100000 )
Write header line for output CSV file that corresponds to stream.
See the comments for run_simulation for a description of the parameters.
"""
function writeheader(stream::IO, simtype::Int64, T::Int64, N_list::Vector{Int64}, N_mu_list::Vector{Float64}, ngens::Int64, 
    burn_in::Float64, slat_reps::Int64=100000 ) 
  param_strings = [
    "# $(string(Dates.today()))",
    "# trials=$(T)",
    "# N_list=$(N_list)",
    "# N_mu_list=$(N_mu_list)",
    "# ngens=$(ngens)",
    "# burn_in=$(burn_in)",
    "# slat_reps=$(slat_reps)",
    ]
  if isdefined(:fract_sample)
    push!(param_strings,"# fract_sample=$(fract_sample)")
  end 
  write(stream, join(param_strings, "\n"), "\n")
  first_heads = ["trial", "N", "N_mu"]
  last_heads = [ "K", "w_homoz", "s_homoz", "s_prob", "e_est_K", "s_est_K" ]
  if simtype == 0
    mid_heads = []
  elseif simtype == 1
    mid_heads = ["cpower"]
  elseif simtype == 2
    mid_heads = ["acer_C", "acer_topsize"]
  end
  line = join(vcat( first_heads, mid_heads, last_heads), ",")
  write(stream, line, "\n")
end

# trial result type that includes Watterson homoz, Slatkin homoz, Slatkin prob, and Ewens K estimate
type trial_result
  N::Int64
  N_mu::Float64
  cpower::Float64
  acer_C::Float64
  acer_topsize::Int64
  K::Int64
  w_homoz::Float64
  s_homoz::Float64
  s_prob::Float64
  ewens_est_K::Float64   #Based on true theta
  s_est_K::Float64       #Based on Slatkin theta est
end

@doc """ function writerow(stream::IOStream, trial::Int64, N::Int64, N_mu::Float64, K::Int64, theta::Float64 )
Write a row to the CSV file.
"""
function writerow(stream::IO, simtype::Int64, trial::Int64, tr::trial_result  )
  first = Any[
    trial,
    tr.N,           # popsize
    tr.N_mu,        # N/mu
  ]
  last = Any[
    tr.K,           # The number of alleles in the sample
    tr.w_homoz,     # Watterson estimated homozygosity
    tr.s_homoz,     # Slatkin estimated homozygosity = 1/(1+theta) where theta = Slatkin est theta
    tr.s_prob,      # Slatkin estimated probability
    tr.ewens_est_K,  # Ewens estimated K based on true theta
    tr.s_est_K,     # Ewens estimated K based Slatkin theta
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
  end
  line = join( vcat( first, mid, last ), "," )
  write(stream, line, "\n")
end

@doc """ function run_trial( N::Int64, N_mu::Float64, ngens::Int64, burn_in::Float64, slat_reps::Int64 )
Run a neutral population evolution and estimate various statitics based on a single generation after burn in.
"""
function run_trial( N::Int64, N_mu::Float64, ngens::Int64, burn_in::Float64=2.0; 
      slat_reps::Int64=10000, cpower::Float64=0.0, acer_C::Float64=0.0, acer_topsize::Int64=10 )
  if cpower==0.0 && acer_C==0.0
    pop = neutral_poplist(N, N_mu/N, ngens, burn_in=burn_in )[ngens]
  elseif cpower != 0.0
    pop = power_conformist_poplist(N, N_mu/N, ngens, burn_in=burn_in, conformist_power=cpower )[ngens]
  elseif acer_C != 0
    pop = acerbi_conformist_poplist(N, N_mu/N, ngens, acer_C, burn_in=burn_in, toplist_size=acer_topsize  )[ngens]
  end
  # If fract_sample is defined in the definition file, take a random subsample
  if isdefined( :fract_sample )
    new_pop = sample_population( pop, Int(ceil( length(pop)*fract_sample )) )
    #println("sampling population with new size: ",Int(ceil( length(pop)*fract_sample)))
  else
    new_pop = pop
  end
  p64 = pop_counts64(new_pop)
  p32 = pop_counts32(new_pop)
  sr = ewens_montecarlo(Int32(slat_reps),p32)
  s_homoz = 1.0/(1.0+sr.theta_estimate)
  w_homoz = watterson_homozygosity(p64)
  ewens_est_K = ewens_K_est( 2*N_mu, N )
  slatkin_est_K = ewens_K_est( sr.theta_estimate, N )
  trial_result(N, N_mu, cpower, acer_C, acer_topsize, length(p64), w_homoz, s_homoz, sr.probability, ewens_est_K, slatkin_est_K )
end

type neutral_trial_type
  N::Int64
  N_mu::Float64
end

type pconform_trial_type
  N::Int64
  N_mu::Float64
  cpower::Float64
end

type aconform_trial_type
  N::Int64
  N_mu::Float64
  acer_C::Float64
  acer_topsize::Int64
end
  
@doc """ function run_simulation(simname::AbstractString, simtype::bool, T::Int64, N_list::Vector{Int64}, N_mu_list::Vector{Float64}, 
    ngens::Int64, burn_in::Float64; slat_reps::Int64==10000 )

This is the function that actually runs the simulation.  
A list of "trials" (arguments to run_trial()) is produced that is fed to pmap so that trials are run in parallel.
Arguments:
simtype:  1 means use Slatkin test
          2 means use Watterson test
          3 means estimate K where K is the number of allele types
T:        the number of trials to run for each setting of N and N_mu
N_list:   a list of population sizes
N_mu_list:  a list of per-generation mutation rates (the per-locus mutation rate is N_mu/N)
ngens:    the number ogenerations to run (not including burn-in)
burn_in:  the number of preliminary generations run to stabilize as a multiple of pop size N
ngens:    the number of generations to run after burn in
slat_reps: the number of monte-carlo reps to use in the Slatkin test
"""

function run_simulation(simname::AbstractString, simtype::Int64, T::Int64, N_list::Vector{Int64}, N_mu_list::Vector{Float64}, 
    ngens::Int64, burn_in::Float64; cpower_list::Vector{Float64}=[0.0], slat_reps::Int64=100000, acer_C_list::Vector{Float64}=[0.0], acer_topsize::Int64=10 )
  stream = open("$(simname).csv", "w")
  if isdefined( :fract_sample )
    multiplier = Int(ceil(1/fract_sample))
    N_list = map(x->multiplier*x, N_list)
    println(" new N_list: ",N_list)
  end
  writeheader(stream,simtype,T,N_list,N_mu_list,ngens,burn_in)
  if simtype == 0   # neutral
    t_list = neutral_trial_type[]
    for N in N_list
      for N_mu in N_mu_list
          push!(t_list,neutral_trial_type(N,N_mu))
      end
    end
  elseif simtype == 1   # power conformity
    t_list = pconform_trial_type[]
    for N in N_list
      for N_mu in N_mu_list
        for cpower in cpower_list
          push!(t_list,pconform_trial_type(N,N_mu,cpower))
        end
      end
    end
  elseif simtype == 2   # Acerbi conformity
    t_list = aconform_trial_type[]
    for N in N_list
      for N_mu in N_mu_list
        for acer_C in acer_C_list
          push!(t_list,aconform_trial_type(N,N_mu,acer_C,acer_topsize))
        end
      end
    end
  end
  trial_list = t_list
  for t in 1:(T-1)
    trial_list = vcat(trial_list,t_list)
  end
  if simtype == 0   # neutral conform
    trials = pmap(tr->run_trial( tr.N, tr.N_mu, ngens, burn_in, slat_reps=slat_reps ), trial_list )
  elseif simtype == 1   # power conform
    trials = pmap(tr->run_trial( tr.N, tr.N_mu, ngens, burn_in, slat_reps=slat_reps, cpower=tr.cpower ), trial_list )
  elseif simtype == 2  # Acerbi conform
    trials = pmap(tr->run_trial( tr.N, tr.N_mu, ngens, burn_in, slat_reps=slat_reps, acer_C=tr.acer_C, acer_topsize=tr.acer_topsize ), trial_list )
  end
  t = 1
  for trial in trials
    writerow(stream, simtype, t, trial )
    flush(stream)
    t += 1
  end
  close(stream)
end

#end #module
