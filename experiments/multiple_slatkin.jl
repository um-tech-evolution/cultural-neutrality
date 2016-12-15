#=
Run multiple Slatkin (or Watterson) tests on successive generations of a 
  simulated infinite-alleles population.
Sample run:
>  julia -p 8
julia> include("multiple_slatkin.jl")   # this file
julia> rtn = setup_neutral_params(50,50,2.0,1)
julia> run_trials_neutral(200,rtn)
=#
@everywhere include("../src/NeutralCulturalEvolution.jl")

using Dates

Btbl = Array{Float64,2}()

@everywhere type neutral_trial_result
  n::Int64
  N::Int64
  N_mu::Float64
  K::Int64
  ngens::Int64
  burn_in::Float64
  sprobs::Vector{Float64}
end

# T is the number of trials
function run_trials_neutral( T::Int64, rtn::neutral_trial_result )
  rtn.sprobs = pmap(x->run_trial_neutral(rtn),zeros(Int64,T))
  write_slat_results(rtn)
end

# T is the number of trials
function run_trials_stewart( T::Int64, stn::neutral_trial_result )
  stn.sprobs = pmap(x->run_trial_stewart(stn),zeros(Int64,T))
  write_slat_results(stn)
end

@everywhere function setup_neutral_params( n::Int64, N::Int64, N_mu::Float64, ngens::Int64=1; burn_in::Float64=2.0 )
  tr = neutral_trial_result(n,N,N_mu,0,ngens,burn_in,Vector{Int64}())
end

@everywhere function setup_stewart_params( n::Int64, N::Int64, K::Int64, ngens::Int64=1; burn_in::Float64=2.0 )
  global Btbl = BT(K,N)
  tr = neutral_trial_result(n,N,0.0,K,ngens,burn_in,Vector{Int64}())
end

icurrent_dir = pwd()
data_string = "../data/"*Dates.format(now(),"mm_dd_yy")*"/"
try mkdir(data_string) catch end   # create today's directory with no error if it already exists
println("data string: ",data_string)

@everywhere function slat_montecarlo( pop::Vector{Int64}, nreps::Int64=100000 )
  ewens_montecarlo( Int32(nreps), pop_counts32(pop) ).probability
end

@doc """ function run_trial_neutral( ntr::neutral_trial_result )
Run a neutral neutral population evolution, then returns the Slatkin probability.
n = sample size
N = popsize
N_mu = mutation rate per population, i. e.,  mu*N.  Thus, mu = N_mu/N.
"""
@everywhere function run_trial_neutral( ntr::neutral_trial_result)
  if ntr.n < ntr.N
    s_m =  slat_montecarlo( sample_population(neutral_poplist(ntr.N, ntr.N_mu, ntr.ngens, burn_in=ntr.burn_in, combine=false )[ntr.ngens],ntr.n))
  else
    s_m = slat_montecarlo(neutral_poplist(ntr.N, ntr.N_mu, ntr.ngens, burn_in=ntr.burn_in, combine=false )[ntr.ngens])
  end
  #filename = "$(data_string)neutral_N:$(N)_N_mu:$(N_mu)_ngens:$(ngens).png"
  #push!(ntr.sprobs, s_m)
  #ntr
  s_m
end

@doc """ function run_trial_stewart( ntr::neutral_trial_result, Btbl::Array{Float64,2} )
Generates a neutral population using the Stewart algorithm, then returns the Slatkin probability.
n = sample size
N = popsize
"""
@everywhere function run_trial_stewart( ntr::neutral_trial_result)
  s_m = slat_montecarlo(Rsample(ntr.N, ntr.K, Btbl ))
end

@doc """ function writeheader(stream::IOStream, t_result::trial_result)
Write all of the trial_result parameters as comments, and then write the header line for Slatkin p-values.
"""
function writeheader(stream::IOStream, t_result::neutral_trial_result)
  write(stream, join([
    "# n=$(t_result.n)",
    "# N=$(t_result.N)",
    "# N_mu=$(t_result.N_mu)",
    "# K=$(t_result.K)",
    "# ngens=$(t_result.ngens)",
    "# burn_in=$(t_result.burn_in)",
    "sprobs"   # the Slatkin p-values from multiple trials
  ], "\n"), "\n")
  #write(stream, line, "\n")
end


@doc """ function write_slat_results( tr::neutral_trial_result,  filename::AbstractString="" )
Write a single-column CSV file where the column is the Slatkin p-values returned by
    multiple runs of run_trial_neutral.
The parameters specified by tr are written as comments preceding the Slatkin p-values.
"""
function write_slat_results( tr::neutral_trial_result,  filename::AbstractString="" )
  filename = "$(data_string)neutral_slatkin_N:$(tr.N)_N_mu:$(tr.N_mu)_ngens:$(tr.ngens).csv"
  stream = open(filename, "w")
  writeheader( stream, tr )
  #write_sprobs( stream, tr.sprobs )
  for i = 1:length(tr.sprobs)
    println(stream, tr.sprobs[i] )
  end
  close(stream)
end

