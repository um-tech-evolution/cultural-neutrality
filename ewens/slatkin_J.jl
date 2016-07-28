export slatkin_sample, slatkin_exact, watterson_sample_gt, watterson_exact_gt, watterson_significance

@doc """ function slatkin_sample( allele_counts::Vector{Int64}, nreps::Int64, Btbl::Array{Float64,2} )
Returns the sum of the probabilities of the allele configurations
   that are less probable that the given configuration.
The result can be used in a two-sided test of the neutral hypothesis.
For example, if the result is less than 0.05, one could reject the null hypothesis at the 0.05 level.
TODO: Write a function that returns both the slatkin and watterson results based on one set of samples.
"""
function slatkin_sample( allele_counts::Vector{Int64}, nreps::Int64, Btbl::Array{Float64,2} )
  N = sum(allele_counts)
  K = length(allele_counts)
  dictc = Dict{Vector{Int64},Int64}()
  dictf = Dict{Vector{Int64},Float64}()
  for j = 1:nreps
    sample = sort(Rsample(N,K,Btbl),rev=true)
    val = get(dictc,sample,0)
    dictc[sample] = val+1
  end
  nrepsf = convert(Float64,nreps)
  for key in keys(dictc)
    dictf[key] = dictc[key]/nrepsf
  end
  pkey = get( dictf, sort(allele_counts,rev=true), 0.0 )
  pkey /= number_orders(allele_counts)
  ssum = 0.0
  for k in keys(dictf)
    if dictf[k]/number_orders(k) <= pkey
      ssum += dictf[k]
    end
  end
  ssum
end

@doc """ function slatkin_sample( allele_counts::Vector{Int64}, nreps::Int64, Btbl::Array{Float64,2} )
Returns the sum of the probabilities of the allele configurations
   that are less probable that the given configuration.
The result can be used in a two-sided test of the neutral hypothesis.
For example, if the result is less than 0.05, one could reject the null hypothesis at the 0.05 level.
TODO: Write a function that returns both the slatkin and watterson results based on one set of samples.
"""
function slatkin_exact( allele_counts::Vector{Int64}, Btbl::Array{Float64,2} )
  N = sum(allele_counts)
  K = length(allele_counts)
  C = cfgs( N, K )[N][K]  # All allele count configs for the given values of N and K
  P = map(c->Pcfg(N,K,c,Btbl),C)   # table of probabilities of configs
  p_alleles = Pcfg(N,K,allele_counts,Btbl)
  ssum = 0.0
  for i = 1:length(C)
    if P[i] <= p_alleles
      ssum += number_orders(C[i])*P[i]
    end
  end
  ssum
end

@doc """ function watterson_sample_gt( allele_counts::Vector{Int64}, nreps::Int64, Btbl::Array{Float64,2} )
Returns the sum of the probabilities of the allele configurations for the given N and K
   whose Watterson homozygosities are greater than or equal to the Watterson homozygosity
   of allele_counts.
The result can be used in a one-sided test of the neutral hypothesis.
For example, if the result is less than 0.05, one could reject the null hypothesis at the 0.05 level.
TODO:  Include multiple options for the direction of the comparison and for one and two sided.
"""
function watterson_sample_gt( allele_counts::Vector{Int64}, nreps::Int64, Btbl::Array{Float64,2} )
  N = sum(allele_counts)
  K = length(allele_counts)
  dictc = Dict{Vector{Int64},Int64}()
  dictf = Dict{Vector{Int64},Float64}()
  for j = 1:nreps
    sample = sort(Rsample(N,K,Btbl),rev=true)
    val = get(dictc,sample,0)
    dictc[sample] = val+1
  end
  nrepsf = convert(Float64,nreps)
  for key in keys(dictc)
    dictf[key] = dictc[key]/nrepsf
  end
  sort!(allele_counts,rev=true)
  pkey = get( dictf, allele_counts, 0.0 )
  pkey /= number_orders(allele_counts)
  wkey =  watterson_homozygosity(allele_counts)
  #println("pkey: ",pkey)
  ssum = 0.0
  for k in keys(dictf)
    if watterson_homozygosity(k) >= wkey
      ssum += dictf[k]
    end
  end
  ssum
end

@doc """ function slatkin_sample( allele_counts::Vector{Int64}, nreps::Int64, Btbl::Array{Float64,2} )
Returns the sum of the probabilities of the allele configurations
   that are less probable that the given configuration.
The result can be used in a two-sided test of the neutral hypothesis.
For example, if the result is less than 0.05, one could reject the null hypothesis at the 0.05 level.
TODO: Write a function that returns both the slatkin and watterson results based on one set of samples.
"""
function watterson_exact_gt( allele_counts::Vector{Int64}, Btbl::Array{Float64,2} )
  N = sum(allele_counts)
  K = length(allele_counts)
  C = cfgs( N, K )[N][K]  # All allele count configs for the given values of N and K
  p_alleles = Pcfg(N,K,allele_counts,Btbl) # probability of given allele_counts 
  w_alleles = watterson_homozygosity( allele_counts )
  P = map(c->Pcfg(N,K,c,Btbl),C)   # table of probabilities of configs
  W = map(watterson_homozygosity,C)   # table of Watterson homozygosities of configs
  ssum = 0.0
  for i = 1:length(C)
    if W[i] >= w_alleles
      ssum += number_orders(C[i])*P[i]
    end
  end
  ssum
end

@doc """ function watterson_significance( N::Int64, K::Int64, Btbl::Array{Float64,2} )
Returns an array of exactly computed significance levels for the Watterson test corresponding to
levels:  0.0, 0.01, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.99, 1.0
Note that the significance level correspoding to 0.0 is the minimum homozygosity, and
the significance level correspoding to 1.0 is the maximum homozygosity.
Agrees fairly closely with Anderson's tables, but not exactly.
"""

function watterson_significance( N::Int64, K::Int64, Btbl::Array{Float64,2} )
  if N > 50 || K > 30
    error("maximum N is 50 and maximum K is 30 in function watterson_significance.")
  end
  C = cfgs( N, K )[N][K]  # All allele count configs for the given values of N and K
  P = map(c->(c,number_orders(c)*Pcfg(N,K,c,Btbl),watterson_homozygosity(c)),C)   # table of (config,probability) tuples 
  #println("P: ",P)
  sort!(P,by=c->c[3])
  #println("P: ",P)
  cummP = deepcopy(P)
  for i in 2:length(cummP)
    cummP[i] = (cummP[i][1],cummP[i][2]+cummP[i-1][2],cummP[i][3])
  end
  levels = [0.0, 0.01, 0.025, 0.05, 0.10, 0.5, 0.90, 0.95, 0.975, 0.99, 1.0]
  ind(s) = findfirst(x->x[2]>s,cummP)
  indices = map(s->ind(s),levels)
  if indices[end] == 0
    indices[end] = length(cummP)
  end
  Float64[ cummP[i][3] for i in indices ]
end

