#include("../src/NeutralCulturalEvolution.jl")


# For power conformist
@doc """ function writeheader_pconform(stream::IOStream, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, ngens::Int64, burn_in::Float64 )
Write header line for output CSV file that corresponds to stream.
See the comments for run_simulation for a description of the parameters.
"""

function writeheader_pconform(stream::IOStream, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, ngens::Int64, burn_in::Float64 )
  write(stream, join([
    "# Power Conformist",
    "# trials=$(T)",
    "# N_list=$(N_list)",
    "# mu_list=$(mu_list)",
    "# ngens=$(ngens)",
    "# burn_in=$(burn_in)",
  ], "\n"), "\n")
  line = join([
    "trial",
    "N",
    "N_mu",
    "cpower",
    "w_theta",
    "s_theta",
    "s_prob",
    "K_est",
    "true_theta",
  ], ",")
  write(stream, line, "\n")
end

# For Acerbi conformist
@doc """ function writeheader_aconform(stream::IOStream, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, ngens::Int64, burn_in::Float64 )
Write header line for output CSV file that corresponds to stream.
See the comments for run_simulation for a description of the parameters.
"""

function writeheader_aconform(stream::IOStream, T::Int64, N_list::Vector{Int64}, mu_list::Vector{Float64}, acer_topsize::Int64, ngens::Int64, burn_in::Float64 )
  write(stream, join([
    "# Acerbi Conformist",
    "# trials=$(T)",
    "# N_list=$(N_list)",
    "# mu_list=$(mu_list)",
    "# acer_topsize=$(acer_topsize)",
    "# ngens=$(ngens)",
    "# burn_in=$(burn_in)",
  ], "\n"), "\n")
  line = join([
    "trial",
    "N",
    "N_mu",
    "acer_C",
    "w_theta",
    "s_theta",
    "s_prob",
    "K_est",
    "true_theta",
  ], ",")
  write(stream, line, "\n")
end

type pconform_trial_result
  N::Int64
  mu::Float64
  cpower::Float64 
  watterson_theta::Float64
  slatkin_theta::Float64
  slatkin_prob::Float64
  K_est::Float64
end

type aconform_trial_result
  N::Int64
  mu::Float64
  acer_topsize::Int64 
  acer_C::Float64 
  watterson_theta::Float64
  slatkin_theta::Float64
  slatkin_prob::Float64
  K_est::Float64
end

# For power conformist
@doc """ function writerow_pconform(stream::IOStream, trial::Int64, trial_result::pconform_trial_result )
Write a row to the CSV file for K estimation.
"""
function writerow_pconform(stream::IOStream, trial::Int64, trial_result::pconform_trial_result )
  line = join(Any[
    trial,
    trial_result.N,
    trial_result.mu,
    trial_result.cpower,
    trial_result.watterson_theta,
    trial_result.slatkin_theta,
    trial_result.slatkin_prob,
    trial_result.K_est,
    2.0*trial_result.mu
  ], ",")
  write(stream, line, "\n")
end

# For Acerbi conformist and anti-conformist
@doc """ function writerow_pconform(stream::IOStream, trial::Int64, trial_result::pconform_trial_result )
Write a row to the CSV file for K estimation.
"""
function writerow_aconform(stream::IOStream, trial::Int64, trial_result::aconform_trial_result )
  line = join(Any[
    trial,
    trial_result.N,
    trial_result.mu,
    trial_result.acer_C,
    trial_result.watterson_theta,
    trial_result.slatkin_theta,
    trial_result.slatkin_prob,
    trial_result.K_est,
    2.0*trial_result.mu
  ], ",")
  write(stream, line, "\n")
end

@doc """ function run_trial_pconform( N::Int64, mu::Float64, ngens::Int64, burn_in::Float64, slat_reps::Int64, cpower::Float64 )  
Run a neutral population evolution and apply K estimation to the result of a single generation.
"""
function run_trial_pconform( N::Int64, mu::Float64, ngens::Int64, burn_in::Float64, slat_reps::Int64, cpower::Float64 )  
  p_counts64 = pop_counts64(power_conformist_poplist(N,mu/N,ngens, burn_in=burn_in, conformist_power=cpower )[ngens])
  p_counts32 = map(x->convert(Int32,x),p_counts64)
  sl = ewens_montecarlo(Int32(slat_reps),p_counts32)
  K_estimate = length(p_counts64)
  w_theta = watterson_theta(p_counts64)
  pconform_trial_result(N, mu, cpower, w_theta, sl.theta_estimate, sl.probability, K_estimate )
end

@doc """ function run_trial_acerbi( N::Int64, mu::Float64, ngens::Int64, burn_in::Float64, C::Float64, toplist_size::Int64 )  
Run a neutral population evolution and apply K estimation to the result of a single generation.
"""
function run_trial_acerbi( N::Int64, mu::Float64, ngens::Int64, burn_in::Float64, slat_reps::Int64, C::Float64, top_list_size::Int64 )  
  p_counts64 = pop_counts64(acerbi_conformist_poplist(N,mu/N,ngens, C, burn_in=burn_in )[ngens])
  p_counts32 = map(x->convert(Int32,x),p_counts64)
  sl = ewens_montecarlo(Int32(slat_reps),p_counts32)
  K_estimate = length(p_counts64)
  w_theta = watterson_theta(p_counts64)
  aconform_trial_result(N, mu, top_list_size, C, w_theta, sl.theta_estimate, sl.probability, K_estimate )
end

#=
N_list = [50,100]
mu_list = [1.0,10.0]
ngens = 1
burn_in = 2.0
cpower = 0.0
C = 0.0
slat_reps = 100000
=#

