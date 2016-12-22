#= Nearly neutral version of simulation
=#

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
#Note:  The location of the shared library file  slatkin.so  must be in LD_LIBRARY_PATH on a Linux system
#@everywhere using NeutralCulturalEvolution  # moved to the end of this file due to what I think is a bug in Julia
include("conformist.jl")
include("../src/NeutralCulturalEvolution.jl")

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
  first_heads = ["trial", "n", "N", "N_mu","K"]
  last_heads = vcat([ "s_prob", 
    #"se_prob", 
    #"we_prob", 
    "s_homoz", "w_homoz"], p_homoz_headers )
  if simtype == 1
    mid_heads = [ "cprob", "acprob", "cpower", "acpower"]
  elseif simtype == 2
    mid_heads = ["acerflg", "bottom", "cprob", "acprob", "topsz", "bottmsz"]
  elseif simtype == 3
    mid_heads = [ "cprob", "acprob", "cpower", "acpower"  ]
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
  cprob::Float64
  acprob::Float64
  cpower::Float64
  acpower::Float64
  topsize::Int64
  bottomsize::Int64
  bottom::Bool
  acerbi_flag::Bool
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
    tr.n,           # sample size
    tr.N,           # popsize
    tr.N_mu,        # N/mu, population mutation rate
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
  elseif simtype == 1 || simtype == 3
    mid = Any[ 
      tr.cprob,
      tr.acprob,
      tr.cpower,
      tr.acpower
    ]
  elseif simtype == 2
    mid = Any[ 
      tr.cprob
      tr.acprob
      tr.bottom
      tr.acerbi_flag
      tr.topsize
      tr.bottomsize
    ]
  end
  line = join( vcat( first, mid, last ), "," )
  write(stream, line, "\n")
end

@doc """ 
Run a population evolution and return various statitics based on a single generation after burn in.
"""
function run_trial( n::Int64, N_mu::Float64, ngens::Int64; psize_mult::Float64=1.0, burn_in::Float64=2.0, slat_reps::Int64=10000, 
      cprob::Float64=0.0, acprob::Float64=0.0, cpower::Float64=0.0, acpower::Float64=0.0, dfe::Function=dfe_neutral,
      acerbi_flag::Bool=true, bottom::Bool=true,topsize::Int64=1,bottomsize::Int64=1)
  if psize_mult > 1.0  # Define popsize N different from sample size n
    N = Int(ceil(psize_mult*n))
    #println("N: ", N)
  else
    N = n
  end
  if simtype == 0
    pop = neutral_poplist(N, N_mu, ngens, burn_in=burn_in, combine=false )[ngens]
  elseif simtype == 1
    pop = power_mixed_conformist_poplist(N, N_mu, ngens, cprob, acprob, conformist_power=cpower, anti_conformist_power=acpower,
      burn_in=burn_in, combine=false )[ngens]
  elseif simtype == 2
    pop = acerbi_mixed_conformist_poplist(N, N_mu, ngens, cprob, acprob, acerbi_flag=true,
      toplist_size=topsize, bottomlist_size=bottomsize, 
      burn_in=burn_in, combine=false  )[ngens]
  elseif simtype == 3
    pop = nearly_neutral_power_mixed_conformist_poplist(N, N_mu, ngens, cprob, acprob, conformist_power=cpower, 
      anti_conformist_power=acpower, dfe = dfe,
      burn_in=burn_in, combine=false )[ngens]
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
  trial_result(N, n, N_mu, length(p64), cprob, acprob, cpower, acpower, topsize, bottomsize, bottom, acerbi_flag, sr.probability, 
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

type pconform_trial_type
  n::Int64
  N_mu::Float64
  cprob::Float64
  acprob::Float64
  cpower::Float64
  acpower::Float64
end

function pconform_trial_type()
  pconform_trial_type(0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

type aconform_trial_type
  n::Int64
  N_mu::Float64
  cprob::Float64
  acprob::Float64
  acerbi_flag::Bool
  bottom::Bool
  topsize::Int64
  bottomsize::Int64
end

function aconform_trial_type()
  aconform_trial_type(0, 0.0, 0.0, 0.0, true, false, 0, 0)
end

type nn_pconform_trial_type
  n::Int64
  N_mu::Float64
  cprob::Float64
  acprob::Float64
  cpower::Float64
  acpower::Float64
  dfe::Function
end

function nn_pconform_trial_type()
  nn_pconform_trial_type(0, 0.0, 0.0, 0.0, 0.0, 0.0, dfe_neutral )
end

function build_global_vars( n_list::Vector{Int64} )
  @everywhere max_n = maximum(n_list)
  @everywhere Btbl = BT(max_n,max_n)
end
  
@doc """ 
This is the function that actually runs the simulation.  
A list of "trials" (arguments to run_trial()) is produced that is fed to pmap so that trials are run in parallel.
Arguments:
simtype:  1 means use power model of conformity (corresponds to length(cpower_list)>0)
          2 means use Acerbi model of conformity 
T:        the number of trials to run for each setting of N and N_mu
n_list:   a list of sample sizes
N_mu_list:  a list of per-generation mutation rates (the per-locus mutation rate is N_mu/N)
ngens:    the number ogenerations to run (not including burn-in)
burn_in:  the number of preliminary generations run to stabilize as a multiple of pop size N
ngens:    the number of generations to run after burn in
slat_reps: the number of monte-carlo reps to use in the Slatkin test
"""
function run_simulation(simname::AbstractString, simtype::Int64, T::Int64, n_list::Vector{Int64}, N_mu_list::Vector{Float64}, 
      ngens::Int64, burn_in::Float64; slat_reps::Int64=100000,  popsize_multiplier::Float64=1.0,
      cprob_list::Vector{Float64}=Float64[], acprob_list::Vector{Float64}=Float64[],
      cpower_list::Vector{Float64}=Float64[], acpower_list::Vector{Float64}=Float64[],
      acer_flag_list::Vector{Bool}=Bool[], bottom_list::Vector{Bool}=Bool[],
      topsize_list::Vector{Int64}=Int64[], bottomsize_list::Vector{Int64}=Int64[],
      dfe::Function=dfe_neutral )
  println("run_simulation")
  for n in n_list
    if n > typemax(ConfigInt)
      error("Reset the type ConfigInt in src/aliases.jl so that n <= typemax(ConfigInt) for all n in n_list.")
    end
  end
  build_global_vars( n_list )
  stream = open("$(simname).csv", "w")
  if popsize_multiplier > 1.0
    N_list = map(x->popsize_multiplier*x, n_list)
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
        for cprob in cprob_list
          for acprob in acprob_list
            for cpower in cpower_list
              for acpower in acpower_list
                push!(t_list,pconform_trial_type(n,N_mu,cprob,acprob,cpower,acpower))
              end
            end
          end
        end
      end
    end
  elseif simtype == 2   # Acerbi conformity
    t_list = aconform_trial_type[]
    for n in n_list
      for N_mu in N_mu_list
        for cprob in cprob_list
          for acprob in acprob_list
            for acer_flag in acer_flag_list
              for bottom in bottom_list
                for topsize in topsize_list
                  for bottomsize in bottomsize_list
                    push!(t_list,aconform_trial_type(n,N_mu,cprob,acprob,acer_flag,bottom,topsize,bottomsize))
                  end
                end
              end
            end
          end
        end
      end
    end
  elseif simtype == 3   # Nearly neutral power conformist
    t_list = nn_pconform_trial_type[]
    for n in n_list
      for N_mu in N_mu_list
        for cprob in cprob_list
          for acprob in acprob_list
            for cpower in cpower_list
              for acpower in acpower_list
                push!(t_list,nn_pconform_trial_type(n,N_mu,cprob,acprob,cpower,acpower,dfe))
              end
            end
          end
        end
      end
    end
  end
  trial_list = t_list
  for t in 1:(T-1)
    trial_list = vcat(trial_list,t_list)
  end
  if simtype == 0   # neutral conform
    trial_results = map(tr->run_trial( tr.n, tr.N_mu, ngens, burn_in=burn_in, psize_mult=popsize_multiplier, slat_reps=slat_reps ), trial_list )
  elseif simtype == 1   # power conform
    trial_results = map(tr->run_trial( tr.n, tr.N_mu, ngens, cprob=tr.cprob, acprob=tr.acprob, cpower=tr.cpower, acpower=tr.acpower, 
        burn_in=burn_in, psize_mult=popsize_multiplier, slat_reps=slat_reps ), trial_list )
  elseif simtype == 2  # Acerbi conform
    trial_results = map(tr->run_trial( tr.n, tr.N_mu, ngens, 
        acerbi_flag=tr.acerbi_flag, bottom=tr.bottom, topsize=tr.topsize, bottomsize=tr.bottomsize,
        burn_in=burn_in, psize_mult=popsize_multiplier, slat_reps=slat_reps ), trial_list )
  elseif simtype == 3  # Nearly neutral power conformist
    trial_results = map(tr->run_trial( tr.n, tr.N_mu, ngens, cprob=tr.cprob, acprob=tr.acprob, cpower=tr.cpower, acpower=tr.acpower, 
        dfe=tr.dfe, burn_in=burn_in, psize_mult=popsize_multiplier, slat_reps=slat_reps ), trial_list )
  end
  trial = 1
  for t_result in trial_results
    writerow(stream, simtype, trial, t_result )
    flush(stream)
    trial += 1
  end
  close(stream)
end
#include("conformist.jl")
#include("../src/NeutralCulturalEvolution.jl")
