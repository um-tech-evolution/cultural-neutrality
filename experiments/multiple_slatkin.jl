#=
Run multiple Slatkin (or Watterson) tests on successive generations of a 
  simulated infinite-alleles population.
Sample run:
>  julia -p 8
julia> include("multiple_slatkin.jl")   # this file
julia> @everywhere rtn = setup_neutral_params(50,50,2.0,1)
julia> run_trials_neutral(200,rtn)
=#
@everywhere include("../src/NeutralCulturalEvolution.jl")

using Dates

@everywhere Btbl = Array{Float64,2}()

@everywhere type neutral_trial_result
  n::Int64
  N::Int64
  N_mu::Float64
  K::Int64
  ngens::Int64
  combine::Bool
  burn_in::Float64
  pvals_no_combine::Vector{Float64}
  pvals_combine::Vector{Float64}
end

# T is the number of trials
function run_trials_neutral( T::Int64, rtn::neutral_trial_result )
  pairs_list = pmap(x->run_trial_neutral(rtn),zeros(Int64,T))
  rtn.pvals_no_combine = [ pairs_list[i][1] for i in 1:T  ]
  rtn.pvals_combine = [ pairs_list[i][2] for i in 1:T  ]
  write_slat_results(rtn)
  rtn
end

# T is the number of trials  # TO DO:  Revise as of 12/21/16
function run_trials_stewart( T::Int64, stn::neutral_trial_result )
  stn.sprobs = pmap(x->run_trial_stewart(stn),zeros(Int64,T))
  write_slat_results(stn)
end

@everywhere function setup_neutral_params( n::Int64, N::Int64, N_mu::Float64, ngens::Int64=1; combine::Bool=false, burn_in::Float64=2.0 )
  tr = neutral_trial_result(n,N,N_mu,0,ngens,combine,burn_in,Vector{Int64}(),Vector{Int64}())
end

@everywhere function setup_stewart_params(  N::Int64, K::Int64, ngens::Int64=1; combine::Bool=false, burn_in::Float64=2.0 )
  global Btbl = BT(K,N)
  tr = neutral_trial_result(0,N,0.0,K,ngens,combine,burn_in,Vector{Int64}(),Vector{Int64}())
end

icurrent_dir = pwd()
data_string = "../data/"*Dates.format(now(),"mm_dd_yy")*"/"
try mkdir(data_string) catch end   # create today's directory with no error if it already exists
println("data string: ",data_string)

@everywhere function slat_montecarlo( pop_counts::Vector{Int64}, nreps::Int64=100000 )
  ewens_montecarlo( nreps, sort(pop_counts,rev=true) ).probability
end

@doc """ function run_trial_neutral( ntr::neutral_trial_result )
Run a neutral neutral population evolution, then returns either the 
    Slatkin probability( if ntr.ngens==1, or the average of the Slatkin probabilities
    of ntr.ngens successive populations
n = sample size
N = popsize
N_mu = mutation rate per population, i. e.,  mu*N.  Thus, mu = N_mu/N.
"""
@everywhere function run_trial_neutral( ntr::neutral_trial_result)
  poplist_no_combine = neutral_poplist(ntr.N, ntr.N_mu, ntr.ngens, burn_in=ntr.burn_in, combine=false )
  #println("p_no_comb:",poplist_no_combine)
  s_m_list = map( x->slat_montecarlo( pop_counts64(sample_population(x,ntr.n))), poplist_no_combine )
  s_m_no_combine = sum(s_m_list)/ntr.ngens
  poplist_combine = combine_pops( poplist_no_combine )
  s_m_combine = slat_montecarlo( pop_counts64(sample_population( poplist_combine, ntr.n )))
  (s_m_no_combine,s_m_combine)
  #=
  if ntr.ngens == 1    # do one generation, return a single p-value
    pop = neutral_poplist(ntr.N, ntr.N_mu, ntr.ngens, burn_in=ntr.burn_in, combine=false )[ntr.ngens]
    s_m =  slat_montecarlo( pop_counts64( sample_population( pop, ntr.n )))
  elseif ntr.ngens > 1 && !ntr.combine  # do ngens generations, average the p-values
    poplist = neutral_poplist(ntr.N, ntr.N_mu, ntr.ngens, burn_in=ntr.burn_in, combine=ntr.combine )[1:ntr.ngens]
    s_m_list = map( x->slat_montecarlo( pop_counts64(sample_population(x,ntr.n))), poplist )
    s_m = sum(s_m_list)/ntr.ngens
  elseif ntr.ngens > 1 && ntr.combine   # do ngens generations, combine the populations, return a single p-value
    poplist = neutral_poplist(ntr.N, ntr.N_mu, ntr.ngens, burn_in=ntr.burn_in, combine=ntr.combine )
    s_m =  slat_montecarlo( pop_counts64( sample_population( pop, ntr.n )))
  end
  s_m
  =#
end

@doc """ function run_trial_stewart( ntr::neutral_trial_result, Btbl::Array{Float64,2} )
Generates a neutral population using the Stewart algorithm, then returns the Slatkin probability.
n = sample size
N = popsize
"""
@everywhere function run_trial_stewart( ntr::neutral_trial_result)
  s_m = slat_montecarlo(sort(Rsample(ntr.N, ntr.K, Btbl ),rev=true))
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
    "p_comb,p_nocomb"   # the Slatkin p-values from combined gens, no combine gens
  ], "\n"), "\n")
  #write(stream, line, "\n")
end


@doc """ function write_slat_results( tr::neutral_trial_result,  filename::AbstractString="" )
Write a single-column CSV file where the column is the Slatkin p-values returned by
    multiple runs of run_trial_neutral.
The parameters specified by tr are written as comments preceding the Slatkin p-values.
"""
function write_slat_results( tr::neutral_trial_result,  filename::AbstractString="" )
  if tr.n == 0 && tr.N_mu == 0.0
    filename = "$(data_string)neutral_slatkin_N:$(tr.N)_N_mu:$(tr.N_mu)_ngens:$(tr.ngens).csv"
  else
    filename = "$(data_string)neutral_inf_alleles_N:$(tr.N)_N_mu:$(tr.N_mu)_ngens:$(tr.ngens).csv"
  end
  stream = open(filename, "w")
  writeheader( stream, tr )
  for i = 1:length(tr.pvals_combine)
    println(stream, tr.pvals_combine[i],",",tr.pvals_no_combine[i] )
  end
  close(stream)
end

