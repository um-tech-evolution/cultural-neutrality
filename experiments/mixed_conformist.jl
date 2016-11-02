#include("../src/NeutralCulturalEvolution.jl")

#=  Copied from conformist.jl on 10/25/16.
Current objective (as of 10/25/16) is to produce a data frame and corresponding CSV file
that can be read into the poweRlaw R package.  So I am deleting the Slatkin and Watterson tests.
=#


# For power conformist
@doc """ function writeheader_pconform(stream::IOStream, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, ngens::Int64, burn_in::Float64 )
Write header line for output CSV file that corresponds to stream.
See the comments for run_simulation for a description of the parameters.
"""

function writeheader_pconform(stream::IOStream, N::Int64, N_mu::Float64, ngens::Int64, cpower::Float64, 
    burn_in::Float64 )
  write(stream, join([
    "# Power Conformist",
    "# N=$(N)",
    "# N_mu=$(N_mu)",
    "# cpower=$(cpower)",
    "# ngens=$(ngens)",
    "# burn_in=$(burn_in)",
  ], "\n"), "\n")
  line = join([
    "pop_counts"
  ], ",")
  write(stream, line, "\n")
end

# For Acerbi conformist
@doc """ function writeheader_aconform(stream::IOStream, T::Int64, N::Int64, mu::Float64, acer_C::Float64, acer_topsize::Int64, ngens::Int64, burn_in::Float64 )
Write header line for output CSV file that corresponds to stream.
See the comments for run_simulation for a description of the parameters.
"""

function writeheader_aconform(stream::IOStream, N::Int64, N_mu::Float64, acer_C::Float64, acer_topsize::Int64, ngens::Int64, burn_in::Float64 )
  write(stream, join([
    "# Acerbi Conformist",
    "# N=$(N)",
    "# N_mu=$(N_mu)",
    "# acer_C=$(acer_C)",
    "# acer_topsize=$(acer_topsize)",
    "# ngens=$(ngens)",
    "# burn_in=$(burn_in)",
  ], "\n"), "\n")
  line = join([
    "pop_counts"
  ], ",")
  write(stream, line, "\n")
end

type neutral_trial_result
  N::Int64
  mu::Float64
  pcounts::Vector{Int64}
end

type pconf_trial_result
  N::Int64
  mu::Float64
  cpower::Float64 
  pcounts::Vector{Int64}
end

type aconf_trial_result
  N::Int64
  mu::Float64
  acer_topsize::Int64 
  acer_C::Float64 
  pcounts::Vector{Int64}
end

type amixed_trial_result
  N::Int64
  mu::Float64
  acer_topsize::Int64 
  conformist_prob::Float64 
  anti_conformist_prob::Float64 
  toplist_size::Int64
  pcounts::Vector{Int64}
end

type pmixed_trial_result
  N::Int64
  mu::Float64
  conformist_prob::Float64 
  anti_conformist_prob::Float64 
  conformist_power::Float64
  anti_conformist_power::Float64
  pcounts::Vector{Int64}
end

# For power conformist
@doc """ function writerow_pconform(stream::IOStream, trial::Int64, trial_result::pconf_trial_result )
Write a row to the CSV file for K estimation.
"""
function writerow_pconform(stream::IOStream, trial_result::pconf_trial_result )
  line = join(Any[
    trial_result.N,
    trial_result.mu,
    trial_result.cpower,
    2.0*trial_result.mu,
    trial_result.pcounts
  ], ",")
  write(stream, line, "\n")
end

# For Acerbi conformist and anti-conformist
@doc """ function writerow_pconform(stream::IOStream, trial::Int64, trial_result::pconf_trial_result )
Write a row to the CSV file for K estimation.
"""
function writerow_aconform(stream::IOStream, trial_result::aconf_trial_result )
  line = join(Any[
    trial,
    trial_result.N,
    trial_result.mu,
    trial_result.acer_C,
    2.0*trial_result.mu,
    trial_result.pcounts
  ], ",")
  write(stream, line, "\n")
end

@doc """ function run_trial_neutral( N::Int64, mu::Float64, ngens::Int64, burn_in::Float64 )  
Run a simple neutral population evolution, computes the poweRlaw parameters, and displays the plot.
n = sample size
N = popsize
"""
#=
function run_trial_neutral( n::Int64,  N::Int64, mu::Float64, ngens::Int64, burn_in::Float64 )  
  if n < N
    p_counts64 = pop_counts64(sample_population(neutral_poplist(N, mu, ngens, burn_in=burn_in )[ngens],n))
  else
    p_counts64 = pop_counts64(neutral_poplist(N, mu, ngens, burn_in=burn_in )[ngens])
  end
  #neutral_trial_result(N, mu, p_counts64 )
  power_law_estimates( p_counts64 )
end
=#
@doc """ function run_trials_simple( n::Int64,  N::Int64, N_mu::Float64, ngens::Int64, burn_in::Float64 )  
Run a simple neutral population evolution combined over ngen generations, 
   then computes the poweRlaw parameters, and displays the plot.
n = sample size
N = popsize
N_mu = mutation rate per population, i. e.,  mu*N.  Thus, mu = N_mu/N.
"""
function run_trials_simple( n::Int64,  N::Int64, N_mu::Float64, ngens::Int64; burn_in::Float64=2.0,
    callR::Bool=true, filename::String="" )  
  if n < N
    p_counts64 = pop_counts64(sample_population(simple_poplist(N, N_mu, ngens, burn_in=burn_in ),n))
  else
    p_counts64 = pop_counts64(simple_poplist(N, N_mu, ngens, burn_in=burn_in ))
  end
  filename = "../data/11_2_16/simple_$(N)_$(N_mu)_$(ngens).png"
  if callR
    power_law_estimates( p_counts64, filename )
  else
    neutral_trial_result(N, N_mu, p_counts64 )
  end
end

@doc """ function run_trials_neutral( n::Int64,  N::Int64, N_mu::Float64, ngens::Int64; burn_in::Float64=2.0, 
    callR::Bool=true )  
Run a neutral neutral population evolution combined over ngen generations, 
   then computes the poweRlaw parameters, and displays the plot.
n = sample size
N = popsize
N_mu = mutation rate per population, i. e.,  mu*N.  Thus, mu = N_mu/N.
"""
function run_trials_neutral( n::Int64,  N::Int64, N_mu::Float64, ngens::Int64; burn_in::Float64=2.0, 
    callR::Bool=true )  
  if n < N
    p_counts64 = pop_counts64(sample_population(neutral_poplist(N, N_mu, ngens, burn_in=burn_in ),n))
  else
    p_counts64 = pop_counts64(neutral_poplist(N, N_mu, ngens, burn_in=burn_in ))
  end
  filename = "../data/11_2_16/neutral_N:$(N)_N_mu:$(N_mu)_ngens:$(ngens).png"
  if callR
    power_law_estimates( p_counts64, filename )
  else
    neutral_trial_result(N, N_mu, p_counts64 )
  end
end

@doc """ function run_trials_pconform( n::Int64, N::Int64, N_mu::Float64, ngens::Int64, cpower::Float64; 
    burn_in::Float64=2.0, callR::Bool=true )
n = sample size
N = popsize
N_mu = mutation rate per population, i. e.,  mu*N.  Thus, mu = N_mu/N.
Run a power conformist population evolution.
"""
function run_trials_pconform( n::Int64, N::Int64, N_mu::Float64, ngens::Int64, cpower::Float64; 
    burn_in::Float64=2.0, callR::Bool=true )
  if n < N
    p_counts64 = pop_counts64(sample_population(power_conformist_poplist(N, N_mu, ngens, 
      conformist_power=cpower, burn_in=burn_in ),n))
  else
    p_counts64 = pop_counts64(power_conformist_poplist(N, N_mu, ngens, conformist_power=cpower, 
      burn_in=burn_in ))
  end
  filename = "../data/11_2_16/pconf_N_$(N)_N_mu_$(N_mu)_ngens_$(ngens)_cpower:$(cpower).png"
  if callR
    power_law_estimates( p_counts64, filename )
  else
    pconf_trial_result(N, N_mu, cpower, p_counts64 )
  end
end

@doc """ function run_trials_acerbi_conformist( n::Int64, N::Int64, N_mu::Float64, ngens::Int64, C::Float64; toplist_size::Int64=1,
    burn_in::Float64=2.0, callR::Bool=true )
n = sample size
N = popsize
N_mu = mutation rate per population, i. e.,  mu*N.  Thus, mu = N_mu/N.
Run an acerbi conformist population evolution.
"""
function run_trials_acerbi_conformist( n::Int64, N::Int64, N_mu::Float64, ngens::Int64, C::Float64; toplist_size::Int64=1,
    burn_in::Float64=2.0, callR::Bool=true )
  if n < N
    p_counts64 = pop_counts64(sample_population(acerbi_conformist_poplist(N, N_mu, ngens, C, 
      toplist_size=toplist_size, burn_in=burn_in ),n))
  else
    p_counts64 = pop_counts64(acerbi_conformist_poplist(N, N_mu, ngens, C, toplist_size=toplist_size, 
      burn_in=burn_in ))
  end
  filename = "../data/11_2_16/acerbi_pos_N_$(N)_N_mu_$(N_mu)_ngens_$(ngens)_C:$(C)_tsize:$(toplist_size).png"
  if callR
    power_law_estimates( p_counts64, filename )
  else
    aconf_trial_result(N, N_mu, toplist_size, C, p_counts64 )
  end
end

@doc """ function run_trials_acerbi( n::Int64, N::Int64, N_mu::Float64, ngens::Int64, C::Float64; toplist_size::Int64=1,
    burn_in::Float64=2.0, callR::Bool=true )
n = sample size
N = popsize
N_mu = mutation rate per population, i. e.,  mu*N.  Thus, mu = N_mu/N.
Run an acerbi conformist population evolution.
"""
function run_trials_acerbi_anti_conformist( n::Int64, N::Int64, N_mu::Float64, ngens::Int64, C::Float64; toplist_size::Int64=1,
    burn_in::Float64=2.0, callR::Bool=true )
  if n < N
    p_counts64 = pop_counts64(sample_population(acerbi_anti_conformist_poplist(N, N_mu, ngens, C, 
      toplist_size=toplist_size, burn_in=burn_in ),n))
  else
    p_counts64 = pop_counts64(acerbi_anti_conformist_poplist(N, N_mu, ngens, C, toplist_size=toplist_size, 
      burn_in=burn_in ))
  end
  filename = "../data/11_2_16/acerbi_neg_N_$(N)_N_mu_$(N_mu)_ngens_$(ngens)_C:$(C)_tsize:$(toplist_size).png"
  if callR
    power_law_estimates( p_counts64, filename )
  else
    aconf_trial_result(N, N_mu, toplist_size, C, p_counts64 )
  end
end

@doc """ function run_trials_acerbi_mixed_conformist( n::Int64, N::Int64, N_mu::Float64, ngens::Int64, 
    conformist_prob::Float64, anti_conformist_prob::Float64; toplist_size::Int64=1,
    burn_in::Float64=2.0, callR::Bool=true )
n = sample size
N = popsize
N_mu = mutation rate per population, i. e.,  mu*N.  Thus, mu = N_mu/N.
Run an acerbi mixed conformist population evolution.
"""
function run_trials_acerbi_mixed_conformist( n::Int64, N::Int64, N_mu::Float64, ngens::Int64, 
    conformist_prob::Float64, anti_conformist_prob::Float64; toplist_size::Int64=1,
    burn_in::Float64=2.0, callR::Bool=true )
  if n < N
    p_counts64 = pop_counts64(sample_population(acerbi_mixed_conformist_poplist(N, N_mu, ngens, 
      conformist_prob, anti_conformist_prob, toplist_size=toplist_size, burn_in=burn_in ),n))
  else
    p_counts64 = pop_counts64(acerbi_mixed_conformist_poplist(N, N_mu, ngens, 
      conformist_prob, anti_conformist_prob, toplist_size=toplist_size, burn_in=burn_in ))
  end
  filename = "../data/11_2_16/acerbi_mixed_N_$(N)_N_mu_$(N_mu)_ngens_$(ngens)_tsize:$(toplist_size)_cprob:$(conformist_prob)_acprob$(anti_conformist_prob).png"
  if callR
    power_law_estimates( p_counts64, filename )
  else
    amixed_trial_result(N, N_mu, toplist_size, conformist_prob, anti_conformist_prob, toplist_size, p_counts64 )
  end
end

@doc """ 
n = sample size
N = popsize
N_mu = mutation rate per population, i. e.,  mu*N.  Thus, mu = N_mu/N.
Run an acerbi mixed conformist population evolution.
"""
function run_trials_power_mixed_conformist( n::Int64, N::Int64, N_mu::Float64, ngens::Int64, 
    conformist_prob::Float64, anti_conformist_prob::Float64; 
    conformist_power::Float64=1.0, anti_conformist_power=1.0,
    burn_in::Float64=2.0, callR::Bool=true )
  if n < N
    p_counts64 = pop_counts64(sample_population(power_mixed_conformist_poplist(N, N_mu, ngens, 
      conformist_prob, anti_conformist_prob, 
      conformist_power=conformist_power, anti_conformist_power=anti_conformist_power,
      burn_in=burn_in ),n))
  else
    p_counts64 = pop_counts64(power_mixed_conformist_poplist(N, N_mu, ngens, 
      conformist_prob, anti_conformist_prob,  
      conformist_power=conformist_power, anti_conformist_power=anti_conformist_power,
      burn_in=burn_in ))
  end
  filename = "../data/11_2_16/power_mixed_N_$(N)_N_mu_$(N_mu)_ngens_$(ngens)_cpower:$(conformist_power)_acpower$(anti_conformist_prob)_cprob:$(conformist_prob)_acprob$(anti_conformist_prob).png"
  if callR
    power_law_estimates( p_counts64, filename )
  else
    pmixed_trial_result(N, N_mu, conformist_prob, anti_conformist_prob, 
      conformist_power, anti_conformist_power, p_counts64 )
  end
end

@doc """ function write_pcounts_pconform( simname::AbstractString, N::Int64, N_mu::Float64, ngens::Int64, 
Write a single-column CSV file (with comments) where the column is the pop counts
    of the population produced by run_trails_pconform.
There is no option for sampling.
"""
function write_pcounts_pconform( simname::AbstractString, N::Int64, N_mu::Float64, ngens::Int64, 
    cpower::Float64; burn_in::Float64=2.0 )
  stream = open("$(simname).csv", "w")
  writeheader_pconform(stream::IOStream, N, N_mu, ngens, cpower, burn_in )
  trial_result = run_trials_pconform( N, N, N_mu, ngens, cpower, burn_in=burn_in, callR=false )
burn_in, 
  #writerow_pconform(stream::IOStream, trial_result::pconf_trial_result )
  write_pcounts( stream, trial_result.pcounts )
  close(stream)
end


function write_pcounts( stream::IOStream, pcounts::Vector{Int64} )
  for i = 1:length(pcounts)
    println(stream, pcounts[i] )
  end
end

#=
n = 50
N = 100
mu = 10.0
ngens = 1
burn_in = 2.0
cpower = 0.0
C = 0.0
=#

