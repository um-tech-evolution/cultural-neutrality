#include("../src/NeutralCulturalEvolution.jl")

#=  Copied from conformist.jl on 10/25/16.
Current objective (as of 10/25/16) is to produce a data frame and corresponding CSV file
that can be read into the poweRlaw R package.  So I am deleting the Slatkin and Watterson tests.
=#

date_string = "../data/11_18_16/"   # make sure that  /home/evotech/cultural_neutrality/data/11_18_16 exists"
println("data string: ",date_string)


@doc """ type trial_result
  A universal trial result for all run_trials functions
  Constructors for specific cases are defined below.
"""
type trial_result
  n::Int64
  N::Int64
  N_mu::Float64
  ngens::Int64
  conformist_prob::Float64 
  anti_conformist_prob::Float64 
  conformist_power::Float64
  anti_conformist_power::Float64
  toplist_size::Int64
  bottomlist_size::Int64
  bottom::Bool
  dfe_str::AbstractString
  burn_in::Float64
  pcounts::Vector{Int64}
end

@doc """ function amixed_trial_result(n::Int64,N::Int64,N_mu::Float64, ngens::Int64, conformist_prob::Float64,
    anti_conformist_prob::Float64, toplist_size::Int64, bottomlist_size::Int64, bottom::Bool,
    burn_in::Float64, pcounts::Vector{Int64} )
Construct a trial_result corresponding to Acerbi mixed
"""
function amixed_trial_result(n::Int64,N::Int64,N_mu::Float64, ngens::Int64, conformist_prob::Float64,
    anti_conformist_prob::Float64, toplist_size::Int64, bottomlist_size::Int64, bottom::Bool,
    burn_in::Float64, pcounts::Vector{Int64} )
  trial_result( n, N, N_mu, ngens, conformist_prob, anti_conformist_prob, 0.0, 0.0, 
      toplist_size, bottomlist_size, bottom, "dfe_none", burn_in, pcounts )
end

@doc """ function pmixed_trial_result(n::Int64,N::Int64,N_mu::Float64, ngens::Int64, conformist_prob::Float64,
    anti_conformist_prob::Float64, conformist_power::Float64, anti_conformist_power::Float64,
    burn_in::Float64, pcounts::Vector{Int64} )
Construct a trial_result corresponding to power mixed
"""
function pmixed_trial_result(n::Int64,N::Int64,N_mu::Float64, ngens::Int64, conformist_prob::Float64,
    anti_conformist_prob::Float64, conformist_power::Float64, anti_conformist_power::Float64,
    burn_in::Float64, pcounts::Vector{Int64} )
  trial_result( n, N, N_mu, ngens, conformist_prob, anti_conformist_prob, conformist_power, anti_conformist_power, 
      0, 0, false, "dfe_none", burn_in, pcounts )
end

function nearly_neutral_trial_result(n::Int64,N::Int64,N_mu::Float64, ngens::Int64, dfe_str::AbstractString,
    burn_in::Float64, pcounts::Vector{Int64} )
  trial_result( n, N, N_mu, ngens, 0.0, 0.0, 0.0, 0.0, 0, 0, false, dfe_str, burn_in, pcounts )
end

@doc """ function writeheader(stream::IOStream, t_result::trial_result)
Write all of the trial_result parameters as comments, and then write the header line for pcounts.
"""
function writeheader(stream::IOStream, t_result::trial_result)
  write(stream, join([
    "# n=$(t_result.n)",
    "# N=$(t_result.N)",
    "# N_mu=$(t_result.N_mu)",
    "# ngens=$(t_result.ngens)",
    "# conformist_prob=$(t_result.conformist_prob)",
    "# anti_conformist_prob=$(t_result.anti_conformist_prob)",
    "# conformist_power=$(t_result.conformist_power)",
    "# anti_conformist_power=$(t_result.anti_conformist_power)",
    "# acer_topsize=$(t_result.toplist_size)",
    "# acer_bottomsize=$(t_result.bottomlist_size)",
    "# bottom=$(t_result.bottom)",
    "# dfe_str=$(t_result.dfe_str)",
    "# burn_in=$(t_result.burn_in)",
  ], "\n"), "\n")
  line = join([
    "pop_counts"
  ], ",")
  write(stream, line, "\n")
end

@doc """ function write_pcounts( filename::AbstractString, tr::trial_result )
Write a single-column CSV file where the column is the pop counts
The parameters specified by tr are written as comments preceding the counts
"""
function write_pcounts( filename::AbstractString, tr::trial_result )
  stream = open(filename, "w")
  writeheader( stream, tr )
  #write_pcounts( stream, tr.pcounts )
  for i = 1:length(tr.pcounts)
    println(stream, tr.pcounts[i] )
  end
  close(stream)
end
@doc """ function run_trials_acerbi_mixed_conformist( n::Int64, N::Int64, N_mu::Float64, ngens::Int64, 
    conformist_prob::Float64, anti_conformist_prob::Float64; toplist_size::Int64=1,
    burn_in::Float64=2.0, callR::Bool=true )
bottom==true means use a bottomlist for anti-conformism
bottom==false means use Acerbi's method of anti-conformism
n = sample size
N = popsize
N_mu = mutation rate per population, i. e.,  mu*N.  Thus, mu = N_mu/N.
Run an acerbi mixed conformist population evolution.
"""
function run_trials_acerbi_mixed_conformist( n::Int64, N::Int64, N_mu::Float64, ngens::Int64, 
    conformist_prob::Float64, anti_conformist_prob::Float64; toplist_size::Int64=1, bottomlist_size::Int64=1,
    bottom::Bool=true, burn_in::Float64=2.0, CSVflag::Bool=true  )
  if n < N
    p_counts64 = pop_counts64(sample_population(acerbi_mixed_conformist_poplist(N, N_mu, ngens, 
      conformist_prob, anti_conformist_prob, toplist_size=toplist_size, bottomlist_size=bottomlist_size, bottom=bottom, burn_in=burn_in ),n))
  else
    p_counts64 = pop_counts64(acerbi_mixed_conformist_poplist(N, N_mu, ngens, 
      conformist_prob, anti_conformist_prob, toplist_size=toplist_size, bottomlist_size=bottomlist_size, bottom=bottom, burn_in=burn_in ))
  end
  if bottom
    filename = "$(date_string)acerbi_mixed_N:$(N)_N_mu:$(N_mu)_ngens:$(ngens)_tsize:$(toplist_size)_bsize:$(bottomlist_size)_cprob:$(conformist_prob)_acprob$(anti_conformist_prob).png"
  else
    filename = "$(date_string)acerbi_mixed_N:$(N)_N_mu:$(N_mu)_ngens:$(ngens)_tsize:$(toplist_size)_cprob:$(conformist_prob)_acprob$(anti_conformist_prob).png"
  end
  println("filename: ",filename)
  if bottom
    filename = "$(date_string)acerbi_mixed_N:$(N)_N_mu:$(N_mu)_ngens:$(ngens)_tsize:$(toplist_size)_bsize:$(bottomlist_size)_cprob:$(conformist_prob)_acprob:$(anti_conformist_prob).png"
  else
    filename = "$(date_string)acerbi_mixed_N:$(N)_N_mu:$(N_mu)_ngens:$(ngens)_tsize:$(toplist_size)_cprob:$(conformist_prob)_acprob:$(anti_conformist_prob).png"
  end
  ple = power_law_estimates( p_counts64, filename )
  atr = amixed_trial_result(n, N, N_mu, ngens, conformist_prob, anti_conformist_prob, 
      toplist_size, bottomlist_size, bottom, burn_in, p_counts64 )
  if CSVflag
    fname = filename[1:(end-4)]*".csv"
    write_pcounts( fname, atr )
  end
  ple
end

@doc """ function run_trials_power_mixed_conformist( n::Int64, N::Int64, N_mu::Float64, ngens::Int64, 
    conformist_prob::Float64, anti_conformist_prob::Float64; 
    conformist_power::Float64=1.0, anti_conformist_power=1.0,
    burn_in::Float64=2.0, CSVflag::Bool=true )
n = sample size
N = popsize
N_mu = mutation rate per population, i. e.,  mu*N.  Thus, mu = N_mu/N.
Run an acerbi mixed conformist population evolution.
"""
function run_trials_power_mixed_conformist( n::Int64, N::Int64, N_mu::Float64, ngens::Int64, 
    conformist_prob::Float64, anti_conformist_prob::Float64; 
    conformist_power::Float64=0.0, anti_conformist_power=0.0,
    burn_in::Float64=2.0, CSVflag::Bool=true )
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
  filename = "$(date_string)power_mixed_N:$(N)_N_mu:$(N_mu)_ngens:$(ngens)_cpower:$(conformist_power)_acpower$(anti_conformist_power)_cprob:$(conformist_prob)_acprob:$(anti_conformist_prob).png"
  ple = power_law_estimates( p_counts64, filename )
  ptr = pmixed_trial_result(n, N, N_mu, ngens, conformist_prob, anti_conformist_prob, 
      conformist_power, anti_conformist_power, burn_in, p_counts64 )
  if CSVflag
    fname = filename[1:(end-4)]*".csv"
    write_pcounts( fname, ptr )
  end
  ple
end

@doc """ 
Run a non-neutral population evolution combined over ngen generations, 
   then computes the poweRlaw parameters, and displays the plot.
n = sample size
N = popsize
N_mu = mutation rate per population, i. e.,  mu*N.  Thus, mu = N_mu/N.
dfe = ditribution of fitness effect function, 
      one of dfe_advantageous, dfe_deleterious, dfe_mixed, dfe_mod
"""
function run_trials_nearly_neutral( n::Int64,  N::Int64, N_mu::Float64, ngens::Int64, dfe::Function; 
    burn_in::Float64=2.0, CSVflag::Bool=true, PNGflag::Bool=true, dfe_mod_fit_inc::Float64=1.0, dfe_mod_modulus::Int64=5 )  
  if dfe == dfe_advantageous
    dfe_str = "adv"
    dfe_funct = dfe_advantageous
  elseif dfe == dfe_deleterious
    dfe_str = "disadv"
    dfe_funct = dfe_deleterious
  elseif dfe == dfe_mixed
    dfe_str = "mixed"
    dfe_funct = dfe_mixed
  elseif dfe == dfe_mod
    if dfe_mod_fit_inc != 1.0
      dfe_str = "mod_fit_inc:$(dfe_mod_fit_inc)modulus:$(dfe_mod_modulus)"
      dfe_funct = x->dfe(x,fit_inc=dfe_mod_fit_inc,modulus=dfe_mod_modulus)
    else
      dfe_str = "mod_neutral"
      dfe_funct = dfe
    end
  else
    dfe_str = "dfe_error"
    dfe_funct = dfe
  end
  if n < N
    p_counts64 = pop_counts64(sample_population(nearly_neutral_poplist(N, N_mu, ngens, dfe_funct, burn_in=burn_in ),n))
  else
    p_counts64 = pop_counts64(nearly_neutral_poplist(N, N_mu, ngens, dfe_funct, burn_in=burn_in ))
  end
  if PNGflag
    filename = "$(date_string)n_neutral_N:$(N)_N_mu:$(N_mu)_ngens:$(ngens)_$(dfe_str).png"
  else  # display plot instead of writing plot file
    filename = ""
  end
  ple = power_law_estimates( p_counts64, filename )
  ntr = nearly_neutral_trial_result(n, N, N_mu, ngens, dfe_str, burn_in, p_counts64 )
  if CSVflag
    fname = filename[1:(end-4)]*".csv"
    write_pcounts( fname, ntr )
  end
  ple
end
#= Set parameters for trial runs
n = N = 250
N_mu = 2.0
ngens = 1000
cprob=0.2
acprob=0.2
tsize=3
bsize=3
cpwr = 0.5
acpwr = -0.5
burn_in = 2.0
dfe = dfe_mod
fit_inc = 1.1
run_trials_power_mixed_conformist(n,N,N_mu,ngens,cprob,acprob,conformist_power=cpwr,anti_conformist_power=acpwr,burn_in=burn_in,CSVflag=true)
run_trials_acerbi_mixed_conformist(n,N,N_mu,ngens,cprob,acprob,toplist_size=tsize,bottomlist_size=bsize,bottom=false,burn_in=burn_in,CSVflag=true)
run_trials_acerbi_mixed_conformist(n,N,N_mu,ngens,cprob,acprob,toplist_size=tsize,bottomlist_size=bsize,bottom=true,burn_in=burn_in,CSVflag=true)
run_trials_nearly_neutral(n,N,N_mu,ngens,dfe,dfe_mod_fit_inc=1.1)
=#

