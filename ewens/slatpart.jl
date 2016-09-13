# Implementation of my algorithm to enumerate the "Slatkin" configurations for a given N and K.
# See the papers "An exact test for neutrality based on Ewens sampling distribution" by Montgomery Slatkin (1994)
#    and "A correction to the exact test for neutrality based on Ewens sampling distribution" (1996).
# N is the sample (or population) size.
# K is the number of distinct alleles
# A configuration is {r_1, r_2, . . . , r_K : r_1 >= r_2 >= . . . >= r_1 >= 1 and n = sum r_i }
# This means that there are r_1 copies of the most frequent allele, r_2 copies of the next most frequent allele, etc.
# For efficiency reasons, I will represent a configuration as a list in reverse order (i. e, smallest first):  
# To avoid confusion, I will denote this as [c_1 <= c_2 <= . . . <= c_K ].  (Thus, c_i = r_{K-i+1}.) 
# ????[r_1, . . . ,r_K],
#   and I will use this notation in the description.
# For each N, there is only 1 configuration with K = 1, namely [N] 
# As explained in the first paper, there are floor(N/2) configurations with K=2, namely
#   [1,N-1], [2,N-2], . . . [floor(N/2),ceil(N/2)].
# Recursively, a configuaration [r_1, r_2, . . . , r_{K-1}, r_K] for N and K, corresponds to a configuration
#   [r_1, r_2, . . . , r_{K-1}] for N-r_K and K-1.
#
# For a bottom up algorithm, a Configig is a list of integers.
# We store configs as a list of lists of lists of configs:
# Given N, K, r_K,  C[N][K][r_K] is a config (which is a list of length K)

export cfgs, ncfgs, rcfgs, flatten, ConfigInt, Config, ConfigList, ConfigConfigList, ConfigConfigConfigList

#=
typealias ConfigInt Int8     # Change to UInt8 if N > 127
#typealias ConfigInt Int128
typealias Config Array{ConfigInt,1}   # A config
typealias ConfigList Array{Array{ConfigInt,1},1}  # A list of configs correponding to N, K 
typealias ConfigConfigList Array{Array{Array{ConfigInt,1},1},1}  # A list of lists of configs corresponding to N
typealias ConfigConfigConfigList Array{Array{Array{Array{ConfigInt,1},1},1},1}  # Corresponds to C
=#

#=
@doc """ function gen_configs()
Generate all configs for up to N = 8.
Entered by hand for debugging purposes.
"""
function gen_configs()
  C1K1 =  Config[ ConfigInt[1] ]    # for N=1, K=1
  C1 = ConfigList[C1K1]    # for N=1
  C = ConfigConfigList[C1]  # all configs
  C2K1 =  Config[ ConfigInt[2] ]   # for N=2, K=1
  C2K2 =  Config[ ConfigInt[1,1] ] # for N=2, K=2
  C2 = ConfigList[C2K1,C2K2]     # for N=2
  push!(C,C2)
  C3K1 = Config[ ConfigInt[3] ]  # for N=3, K=1
  C3K2 = Config[ ConfigInt[1,2] ]  # for N=3, K=2
  C3K3 = Config[ ConfigInt[1,1,1] ]  # for N=3, K=3
  C3 = ConfigList[C3K1,C3K2,C3K3]     # for N=3
  push!(C,C3)
  C4K1 = Config[ ConfigInt[4] ]  # for N=4, K=1
  C4K2 = Config[ ConfigInt[2,2], ConfigInt[1,3] ]  # for N=4, K=2  
  C4K3 = Config[ ConfigInt[1,1,2] ]  # for N=4, K=3
  C4K4 = Config[ ConfigInt[1,1,1,1] ]  # for N=4, K=4
  C4 = ConfigList[C4K1,C4K2,C4K3,C4K4]     # for N=4
  push!(C,C4)
  C5K1 = Config[ ConfigInt[5] ]  # for N=5, K=1
  C5K2 = Config[ ConfigInt[2,3], ConfigInt[1,4]  ]  # for N=5, K=2  ?
  C5K3 = Config[ ConfigInt[1,2,2], ConfigInt[1,1,3] ]  # for N=5, K=3
  C5K4 = Config[ ConfigInt[1,1,1,2] ]  # for N=5, K=4
  C5K5 = Config[ ConfigInt[1,1,1,1,1] ]  # for N=5, K=5
  C5 = ConfigList[C5K1,C5K2,C5K3,C5K4,C5K5]     # for N=5
  push!(C,C5)
  C6K1 = Config[ ConfigInt[6] ]  # for N=6, K=1
  C6K2 = Config[ ConfigInt[3,3], ConfigInt[2,4], ConfigInt[1,5] ]  # for N=6, K=2  ?
  C6K3 = Config[ ConfigInt[2,2,2], ConfigInt[1,2,3], ConfigInt[1,1,4] ]  # for N=6, K=3
  C6K4 = Config[ ConfigInt[1,1,2,2], ConfigInt[1,1,1,3] ]  # for N=6, K=4
  C6K5 = Config[ ConfigInt[1,1,1,1,2] ]  # for N=6, K=5
  C6K6 = Config[ ConfigInt[1,1,1,1,1,1] ]  # for N=6, K=6
  C6 = ConfigList[C6K1,C6K2,C6K3,C6K4,C6K5,C6K6]     # for N=6
  push!(C,C6)
  C7K1 = Config[ ConfigInt[7] ]  # for N=7, K=1
  C7K2 = Config[ ConfigInt[3,4], ConfigInt[2,5], ConfigInt[1,6] ]  # for N=7, K=2  ?
  C7K3 = Config[ ConfigInt[2,2,3], ConfigInt[1,3,3], ConfigInt[1,2,4], ConfigInt[1,1,5] ]  # for N=7, K=3
  C7K4 = Config[ ConfigInt[1,2,2,2],  ConfigInt[1,1,2,3],  ConfigInt[1,1,1,4] ]  # for N=7, K=4
  C7K5 = Config[ ConfigInt[1,1,1,2,2], ConfigInt[1,1,1,1,3] ]  # for N=7, K=5
  C7K6 = Config[ ConfigInt[1,1,1,1,1,2] ]  # for N=7, K=6
  C7K7 = Config[ ConfigInt[1,1,1,1,1,1,1] ]  # for N=7, K=7
  C7 = ConfigList[C7K1,C7K2,C7K3,C7K4,C7K5,C7K6,C7K7]     # for N=7
  push!(C,C7)
  C8K1 = Config[ ConfigInt[8] ]  # for N=8, K=1
  C8K2 = Config[ ConfigInt[4,4], ConfigInt[3,5], ConfigInt[2,6], ConfigInt[1,7] ]  # for N=8, K=2  ?
  C8K3 = Config[ ConfigInt[2,3,3], ConfigInt[2,2,4], ConfigInt[1,3,4], ConfigInt[1,2,5], ConfigInt[1,1,6] ]  # for N=8, K=3
  C8K4 = Config[ ConfigInt[2,2,2,2],  ConfigInt[1,2,2,3], ConfigInt[1,1,3,3], ConfigInt[1,1,2,4], ConfigInt[1,1,1,5] ]  # for N=8, K=4
  C8K5 = Config[ ConfigInt[1,1,2,2,2], ConfigInt[1,1,1,2,3], ConfigInt[1,1,1,1,4] ]  # for N=8, K=5
  C8K6 = Config[ ConfigInt[1,1,1,1,2,2], ConfigInt[1,1,1,1,1,3] ]  # for N=8, K=6
  C8K7 = Config[ ConfigInt[1,1,1,1,1,1,2] ]  # for N=8, K=7
  C8K8 = Config[ ConfigInt[1,1,1,1,1,1,1,1] ]  # for N=8, K=8
  C8 = ConfigList[C8K1,C8K2,C8K3,C8K4,C8K5,C8K6,C8K7,C8K8]     # for N=8
  push!(C,C8)
  return C
end
=#

@doc """ function ncfgs( N::Int64, K::Int64 )
Call rcfgs to Generate the configs for N and K recursively
"""
function ncfgs( N::Int64, K::Int64 )
  if K > N || K < 1
    error("illegal value of K in function ncfgs()")
  end
  result= Vector{ConfigInt}[ConfigInt[0 for i = 1:K]]
  rcfgs( N, K, result )   # The last config is duplicated
  pop!(result)
  result
end

@doc """ 
Recursively generate all configs for a given n and k.
"""
function rcfgs( n::Int64, k::Int64, result::Array{Vector{ConfigInt},1} )
  #println("cfgs: n:",n,"  k:",k)
  K = length(result[1])
  maxj = (k == K) ? n-k+1 : min(result[end][k+1],n-k+1)
  #println("n:",n,"  k:",k,"  minj: ",cld(n,k),"  maxj: ",maxj)
  for j = cld(n,k):maxj
    result[end][k] = j
    #println("n:",n,"  k:",k,"  j: ",j,"  res:",transpose(result[end]))
    if k == 1
      nresult = deepcopy(result[end])
      #println("nresult: ",transpose(nresult))
      push!(result,nresult)
    end
    if k > 1
      rcfgs(n-j,k-1,result)
    end
  end
end

function rcfs( n::Int64, k::Int64, result::Array{Vector{ConfigInt},1} )
  #println("cfgs: n:",n,"  k:",k)
  K = length(result[1])
  maxj = (k == K) ? n-k+1 : min(result[end][k+1],n-k+1)
  #println("n:",n,"  k:",k,"  minj: ",cld(n,k),"  maxj: ",maxj)
  for j = cld(n,k):maxj
    result[end][k] = j
    #println("n:",n,"  k:",k,"  j: ",j,"  res:",transpose(result[end]))
    if k == 2
      tmplist = [ ConfigInt[i,n-i] for i = ConfigInt(div(n,2)):ConfigInt(-1):ConfigInt(1)]
      result = vcat(result[1:(end-1)],tmplist)
      println("k=2: ",result)
      return result
    end
    if k == 1 
      nresult = deepcopy(result[end])
      println("k=1:  result: ",transpose(result))
      push!(result,nresult)
      return result
    end
    if k > 2
      rcfs(n-j,k-1,result)
    end
  end
end


@doc """ function cfgs( N::Int64, K::Int64 )
Generate the complete table of configs for up to N and K using a bottom-up algorithm.
"""
function cfgs( N::Int64, K::Int64 )
  C1K1 =  Config[ ConfigInt[1] ]    # for N=1, K=1
  C1 = ConfigList[C1K1]    # for N=1
  C = ConfigConfigList[C1]  # all configs
  for n = 2:N
    nresult = ConfigList[ Config[ ConfigInt[n] ] ]
    for k = 2:min(n,K)
      kresult = Config[]
      for j = cld(n,k):(n-k+1)
        fff = C[n-j][k-1]
        ggg = extnd( fff, j)
        for c in ggg
          push!(kresult,c)
        end
      end
      #println("n:",n,"  k:",k,"  kresult: ",kresult)
      push!(nresult,kresult)
    end
    #println("n:",n,"  nresult: ",nresult)
    push!(C,nresult)
  end
  C
end

#= Replace by rcfgs() above which is faster and more memory efficient
@doc """ function rcfgs( N::Int64, K::Int64 )
Recursively generate all configs for a given N and K.
"""
function rcfgs( N::Int64, K::Int64 )
  #println("cfgs: N:",N,"  K:",K)
  if K > N || K < 1
    error("illegal value of K in function cfgs()")
  end
  if K == 1
    return Config[ConfigInt[N]]
  end
  result = Config[]
  for j = cld(N,K):(N-K+1)
    fff = rcfgs(N-j,K-1)
    ggg = extnd( fff, j)
    for c in ggg
      push!(result,c)
    end
  end
  result
end
=#

@doc """ function extnd( lst::ConfigList, elt::Int64 )
"""
function extnd( allele_counts::ConfigList, elt::Int64 )
  result = Config[]
  for r in allele_counts
    if r[end] <= elt
      s = deepcopy(r)
      push!(s,elt)
      push!(result,s)
    end
  end
  result
end

@doc """ function flatten( lst::Array{Array{Array{Int64,1},1},1} )
Flattens a list of list of allele count configs into a list of allele count configs.
Example:  
julia> C = cfgs(4,4);  flatten(C[4])
5-element Array{Array{Int64,1},1}:
 [4]
 [2,2]
 [1,3]
 [1,1,2]
 [1,1,1,1]
This is a list of all allele count configurations of length 4.
"""
function flatten( lst::Array{Array{Array{Int64,1},1},1})
  result = Array{Int64,1}[]
  for p in lst
    for q in p
      push!(result,q)
    end
  end
  result
end
