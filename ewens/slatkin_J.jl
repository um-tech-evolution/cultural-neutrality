export Rsamples, slat_sample, slatkin_sample, slatkin_exact, watterson_sample, watterson_exact, watterson_significance

# isless is aliased to < by default
#typealias islessequal <=
#typealias isgreater >
#typealias isgreaterequal >=
const islessequal = <=
const isgreater = >
const isgreaterequal = >=
#typealias configs_type Array{Array{Array{Array{Int8,1},1},1},1}
#typealias ConfigList Array{Array{Int8,1},1}

immutable rsamples_result
  dict::Dict{Config,Float64}
  ave_whomozyg::Float64
end

@doc """ Rsamples( N::Int64, K::Int64, nreps::Int64, Btbl::Array{Float64,2} )
Returns a pair where the first element of the pair is a dictionary of allele count samples
  computed by the Rsample function (which implements Stewart's sampling algorithm),
  and the second element of the pair is the average Watterson homozygosity of the samples.
"""
function Rsamples( N::Int64, K::Int64, nreps::Int64, Btbl::Array{Float64,2} )
  dictc = Dict{Config,Int64}()
  dictf = Dict{Config,Float64}()
  for j = 1:nreps
    sample = sort(Rsample(N,K,Btbl),rev=true)
    val = get(dictc,sample,0)
    dictc[sample] = val+1
  end
  nrepsf = convert(Float64,nreps)
  wsum = 0.0   # sum Watterson homozygosities
  for cfg in keys(dictc)
    wsum += watterson_homozygosity(cfg)
    dictf[cfg] = dictc[cfg]/nrepsf
  end
  average_watterson_homozygosity = wsum/nrepsf
  rsamples_result(dictf, average_watterson_homozygosity)
end

@doc """ function slat_sample( allele_counts::Config, nreps::Int64, C::ConfigList, Btbl::Array{Float64,2} )
Correctly approximates the Slatkin exact test probability of allele_counts
by taking a random unbiased sample of the ConfigList of all configs.
But impractical because it requires the list of all configs for the given sample size n and config length k.
"""
function slat_sample( allele_counts::Config, nreps::Int64, Btbl::Array{Float64,2}, C::ConfigList )
  if nreps > 0.9*length(C)
    error("nreps too large in r_slat")
  end
  Cind_set = Set{Int64}()
  #Csublist = zeros( Config, nreps )
  pkey = Pcfg( allele_counts, Btbl )
  total_sum = 0.0
  ssum = 0.0
  for i = 1:nreps
    r = rand(1:length(C))
    while r in Cind_set
      r = rand(1:length(C))
    end
    cfg = C[r]
    pc = Pcfg( cfg, Btbl )
    total_sum += pc*number_ordersF(cfg)
    if pc < pkey
      ssum += number_ordersF(cfg)*pc
    end
  end
  ssum/total_sum
end

@doc """ function slatkin_sample( allele_counts::Config, samples::rsamples_result )
Returns the sum of the probabilities of the allele configurations
   that are less probable that the given configuration.
The result can be used in a two-sided test of the neutral null hypothesis.
For example, if the result is less than 0.025 or greater than 0.975, one could reject the null hypothesis 
at the 0.05 level.
TODO:  Add the Watterson probability which is included in Slatkin's code.
"""
function slatkin_sample( allele_counts::Vector{Int64}, nreps::Int64, Btbl::Array{Float64,2} )
  N = sum(allele_counts)
  K = length(allele_counts)
  w_allele_counts =  watterson_homozygosity(allele_counts)
  pkey = 1.0/prod(allele_counts)
  e_count = 0
  f_count = 0
  for j = 1:nreps
    sample = Rsample(N,K,Btbl)
    if 1.0/prod(sample) <= pkey
      e_count += 1
    end
    if watterson_homozygosity(sample) <= w_allele_counts
      f_count += 1
    end
  end
  nrepsf = convert(Float64,nreps)
  (e_count/nrepsf, f_count/nrepsf)
end

@doc """ function slatkin_exact( allele_counts::Config, Btbl::Array{Float64,2}, C::ConfigList=ConfigList() )
Returns the sum of the probabilities of the allele configurations
   that are less probable or equal that the given configuration.
The result can be used in a two-sided test of the neutral null hypothesis.
For example, if the result is less than 0.025 or greater than 0.975, one could reject the null hypothesis 
at the 0.05 level.
TODO: Write a function that returns both the slatkin and watterson results based on one set of samples.
"""
function slatkin_exact( allele_counts::Config, Btbl::Array{Float64,2}, C::ConfigList=ConfigList() )
  N = convert(Int64,sum(allele_counts))
  K = length(allele_counts)
  if length(C) == 0
    C = ncfgs( N, K )    # All allele count configs for the given values of N and K
  end
  P = map(c->Pcfg(c,Btbl),C)   # table of probabilities of configs
  p_alleles = Pcfg(N,K,allele_counts,Btbl)
  ssum = 0.0
  for i = 1:length(C)
    if P[i] <= p_alleles
      ssum += number_ordersF(C[i])*P[i]
    end
  end
  ssum
end

@doc """ function watterson_sample( allele_counts::Config, samples::rsamples_result;
    lt3::Function=islessequal )
Returns the sum of the probabilities of the allele configurations cfg for the given N and K
  such that watterson_homozygosity(cfg) <= watterson_homozygosity(allele_counts).  (dfeault comparison) 
The direction of comparison can be changed by specifying the comparison function lte.
The result can be used in a one or two-sided test of the neutral hypothesis.
For example, if the result is greater than 0.95, one could reject the null hypothesis for the one-sided
   alternative that the homozygosity is too large to be neutral at the 0.05 level.
This function is inaccurate for large N because configurations with low probability are not sufficiently sampled.
"""
function watterson_sample( allele_counts::Config;
    lte::Function=islessequal )
  N = sum(allele_counts)
  K = length(allele_counts)
  w_allele_counts =  watterson_homozygosity(allele_counts)
  pkey = get( samples.dict, sort(allele_counts,rev=true), 0.0 )
  pkey /= number_ordersF(allele_counts)
  w_allele_counts =  watterson_homozygosity(allele_counts)
  #println("pkey: ",pkey)
  ssum = 0.0
  for cfg in keys(samples.dict)
    if lte( watterson_homozygosity(cfg), w_allele_counts )
      ssum += samples.dict[cfg]
    end
  end
  ssum
end

@doc """ function watterson_exact( allele_counts::Config, Btbl::Array{Float64,2}, C::ConfigList=ConfigList(); lte::Function=islessequal )
Returns the sum of the probabilities of the allele configurations cfg for the given N and K
  such that watterson_homozygosity(cfg) <= watterson_homozygosity(allele_counts).  (dfeault comparison) 
The direction of comparison can be changed by specifying the comparison function lte.
The result can be used in a one or two-sided test of the neutral hypothesis.
For example, if the result is greater than 0.95, one could reject the null hypothesis for the one-sided
   alternative that the homozygosity is too large to be neutral at the 0.05 level.
TODO: Write a function that returns both the slatkin and watterson results based on one set of samples.
"""
function watterson_exact( allele_counts::Config, Btbl::Array{Float64,2}, C::ConfigList=ConfigList(); lte::Function=islessequal )
  N = sum(allele_counts)
  K = length(allele_counts)
  if size(C)[1] == 0
    C = ncfgs( N, K )    # All allele count configs for the given values of N and K
  end
  #C = cfgs( N, K )[N][K]  # All allele count configs for the given values of N and K
  p_alleles = Pcfg(N,K,allele_counts,Btbl) # probability of given allele_counts 
  w_alleles = watterson_homozygosity( allele_counts )
  P = map(c->Pcfg(N,K,c,Btbl),C)   # table of probabilities of configs
  W = map(watterson_homozygosity,C)   # table of Watterson homozygosities of configs
  ssum = 0.0
  for i = 1:length(C)
    if lte( W[i], w_alleles )
      ssum += number_ordersF(C[i])*P[i]
    end
  end
  ssum
end

@doc """ function watterson_significance( N::Int64, K::Int64, Btbl::Array{Float64,2} )
Returns an array of exactly computed significance levels for the Watterson test corresponding to
levels:  0.0, 0.01, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.99, 1.0
Note that the significance level correspoding to 0.0 is the minimum homozygosity, and
the significance level correspoding to 1.0 is the maximum homozygosity.
If the p parameter is different from 2, then uses p-significance rather than Watterson significance.
Agrees fairly closely with Anderson's tables, but not exactly.
"""

function watterson_significance( N::Int64, K::Int64, Btbl::Array{Float64,2}; nreps::Int64=0, 
    C::ConfigList=ConfigList(), p::Float64=2.0 )
  #=
  if nreps == 0 && (N > 50 || K > 30)
    error("maximum N is 50 and maximum K is 30 in function watterson_significance using exact computation (nreps==0).")
  end
  =#
  if nreps === 0 && size(C)[1] == 0
    C = cfgs( N, K )[N][K]  # All allele count configs for the given values of N and K
  elseif nreps > 0
    samples = Rsamples(N,K,nreps,Btbl)
    C = [ cfg for cfg in keys(samples.dict) ]
  end
  if p == 2.0
    P = map(c->(c,number_ordersF(c)*Pcfg(N,K,c,Btbl),watterson_homozygosity(c)),C)   # table of (config,probability) tuples 
  else
    P = map(c->(c,number_ordersF(c)*Pcfg(N,K,c,Btbl),p_homozygosity(c,p)),C)   # table of (config,probability) tuples 
  end
  #println("P: ",P)
  sort!(P,by=c->c[3])
  #println("P: ",P)
  cummP = deepcopy(P)
  for i in 2:length(cummP)
    cummP[i] = (cummP[i][1],cummP[i][2]+cummP[i-1][2],cummP[i][3])
  end
  if cummP[end][2] < 1.0-10.0*eps()
    multiplier = 1.0/cummP[end][2]
    println("normalizing multiplier: ",multiplier)
    for i in 1:length(cummP)
      cummP[i] = (cummP[i][1],multiplier*cummP[i][2],cummP[i][3])
    end
  end
  levels = [0.0, 0.01, 0.025, 0.05, 0.10, 0.5, 0.90, 0.95, 0.975, 0.99, 1.0]
  s_points = Any[ 0 for i = 1:length(levels) ]
  s_points[1] = cummP[1][3]      # min homozygosity
  s_points[end] = cummP[end][3]  # max homozygosity
	half_index = findfirst(x->x>=0.5,levels)
  for i = 2:half_index
		j = findlast(x->x[2]<levels[i],cummP)
		if j == 0
      s_points[i] = "-"
		else 
      s_points[i] = cummP[j][3]
		end
  end
  for i = half_index+1:length(levels)-1
		j = findfirst(x->x[2]>levels[i],cummP)
		#println("i: ",i,"  j: ",j)
		if j == 0 || j == length(cummP)
      s_points[i] = "-"
		else 
      s_points[i] = cummP[j][3]
		end
  end
  (s_points, cummP )
end
