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
# For a bottom up algorithm, a Config is a list of integers.
# We store configs as a list of lists of lists of configs:
# Given N, K, r_K,  C[N][K][r_K] is a config (which is a list of length K)

export cfgs, rcfgs, flatten

typealias MyInt Int64
#typealias MyInt Int128
typealias Conf Array{MyInt,1}   # A config
typealias ConfList Array{Array{MyInt,1},1}  # A list of configs correponding to N, K 
typealias ConfConfList Array{Array{Array{MyInt,1},1},1}  # A list of lists of configs corresponding to N
typealias ConfConfConfList Array{Array{Array{Array{MyInt,1},1},1},1}  # Corresponds to C

#=
@doc """ function gen_configs()
Generate all configs for up to N = 8.
Entered by hand for debugging purposes.
"""
function gen_configs()
  C1K1 =  Conf[ MyInt[1] ]    # for N=1, K=1
  C1 = ConfList[C1K1]    # for N=1
  C = ConfConfList[C1]  # all configs
  C2K1 =  Conf[ MyInt[2] ]   # for N=2, K=1
  C2K2 =  Conf[ MyInt[1,1] ] # for N=2, K=2
  C2 = ConfList[C2K1,C2K2]     # for N=2
  push!(C,C2)
  C3K1 = Conf[ MyInt[3] ]  # for N=3, K=1
  C3K2 = Conf[ MyInt[1,2] ]  # for N=3, K=2
  C3K3 = Conf[ MyInt[1,1,1] ]  # for N=3, K=3
  C3 = ConfList[C3K1,C3K2,C3K3]     # for N=3
  push!(C,C3)
  C4K1 = Conf[ MyInt[4] ]  # for N=4, K=1
  C4K2 = Conf[ MyInt[2,2], MyInt[1,3] ]  # for N=4, K=2  
  C4K3 = Conf[ MyInt[1,1,2] ]  # for N=4, K=3
  C4K4 = Conf[ MyInt[1,1,1,1] ]  # for N=4, K=4
  C4 = ConfList[C4K1,C4K2,C4K3,C4K4]     # for N=4
  push!(C,C4)
  C5K1 = Conf[ MyInt[5] ]  # for N=5, K=1
  C5K2 = Conf[ MyInt[2,3], MyInt[1,4]  ]  # for N=5, K=2  ?
  C5K3 = Conf[ MyInt[1,2,2], MyInt[1,1,3] ]  # for N=5, K=3
  C5K4 = Conf[ MyInt[1,1,1,2] ]  # for N=5, K=4
  C5K5 = Conf[ MyInt[1,1,1,1,1] ]  # for N=5, K=5
  C5 = ConfList[C5K1,C5K2,C5K3,C5K4,C5K5]     # for N=5
  push!(C,C5)
  C6K1 = Conf[ MyInt[6] ]  # for N=6, K=1
  C6K2 = Conf[ MyInt[3,3], MyInt[2,4], MyInt[1,5] ]  # for N=6, K=2  ?
  C6K3 = Conf[ MyInt[2,2,2], MyInt[1,2,3], MyInt[1,1,4] ]  # for N=6, K=3
  C6K4 = Conf[ MyInt[1,1,2,2], MyInt[1,1,1,3] ]  # for N=6, K=4
  C6K5 = Conf[ MyInt[1,1,1,1,2] ]  # for N=6, K=5
  C6K6 = Conf[ MyInt[1,1,1,1,1,1] ]  # for N=6, K=6
  C6 = ConfList[C6K1,C6K2,C6K3,C6K4,C6K5,C6K6]     # for N=6
  push!(C,C6)
  C7K1 = Conf[ MyInt[7] ]  # for N=7, K=1
  C7K2 = Conf[ MyInt[3,4], MyInt[2,5], MyInt[1,6] ]  # for N=7, K=2  ?
  C7K3 = Conf[ MyInt[2,2,3], MyInt[1,3,3], MyInt[1,2,4], MyInt[1,1,5] ]  # for N=7, K=3
  C7K4 = Conf[ MyInt[1,2,2,2],  MyInt[1,1,2,3],  MyInt[1,1,1,4] ]  # for N=7, K=4
  C7K5 = Conf[ MyInt[1,1,1,2,2], MyInt[1,1,1,1,3] ]  # for N=7, K=5
  C7K6 = Conf[ MyInt[1,1,1,1,1,2] ]  # for N=7, K=6
  C7K7 = Conf[ MyInt[1,1,1,1,1,1,1] ]  # for N=7, K=7
  C7 = ConfList[C7K1,C7K2,C7K3,C7K4,C7K5,C7K6,C7K7]     # for N=7
  push!(C,C7)
  C8K1 = Conf[ MyInt[8] ]  # for N=8, K=1
  C8K2 = Conf[ MyInt[4,4], MyInt[3,5], MyInt[2,6], MyInt[1,7] ]  # for N=8, K=2  ?
  C8K3 = Conf[ MyInt[2,3,3], MyInt[2,2,4], MyInt[1,3,4], MyInt[1,2,5], MyInt[1,1,6] ]  # for N=8, K=3
  C8K4 = Conf[ MyInt[2,2,2,2],  MyInt[1,2,2,3], MyInt[1,1,3,3], MyInt[1,1,2,4], MyInt[1,1,1,5] ]  # for N=8, K=4
  C8K5 = Conf[ MyInt[1,1,2,2,2], MyInt[1,1,1,2,3], MyInt[1,1,1,1,4] ]  # for N=8, K=5
  C8K6 = Conf[ MyInt[1,1,1,1,2,2], MyInt[1,1,1,1,1,3] ]  # for N=8, K=6
  C8K7 = Conf[ MyInt[1,1,1,1,1,1,2] ]  # for N=8, K=7
  C8K8 = Conf[ MyInt[1,1,1,1,1,1,1,1] ]  # for N=8, K=8
  C8 = ConfList[C8K1,C8K2,C8K3,C8K4,C8K5,C8K6,C8K7,C8K8]     # for N=8
  push!(C,C8)
  return C
end
=#

@doc """ function cfgs( N::Int64, K::Int64 )
Generate the complete table of configs for up to N and K using a bottom-up algorithm.
"""
function cfgs( N::Int64, K::Int64 )
  C1K1 =  Conf[ MyInt[1] ]    # for N=1, K=1
  C1 = ConfList[C1K1]    # for N=1
  C = ConfConfList[C1]  # all configs
  for n = 2:N
    nresult = ConfList[ Conf[ MyInt[n] ] ]
    for k = 2:min(n,K)
      kresult = Conf[]
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

@doc """ function rcfgs( N::Int64, K::Int64 )
Recursively generate all configs for a given N and K.
"""
function rcfgs( N::Int64, K::Int64 )
  #println("cfgs: N:",N,"  K:",K)
  if K > N || K < 1
    error("illegal value of K in function cfgs()")
  end
  if K == 1
    return Conf[MyInt[N]]
  end
  result = Conf[]
  for j = cld(N,K):(N-K+1)
    fff = rcfgs(N-j,K-1)
    ggg = extnd( fff, j)
    for c in ggg
      push!(result,c)
    end
  end
  result
end

@doc """ function extnd( lst::ConfList, elt::Int64 )
"""
function extnd( allele_counts::ConfList, elt::Int64 )
  result = Conf[]
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

# TODO: delete or move to test
function test_all(N::Int64,C::ConfConfConfList)
  for n = 2:N
    for k = 1:n
      cfg = rcfgs(n,k)
      println("n: ",n,"  k: ",k,"  cfg: ",cfg)
      if cfg != C[n][k]
        println("n: ",n,"  k: ",k,"  cfg: ",cfg,"  C: ",C[n][k])
      end
    end
  end
end

#= TODO:  delete these three functions
function Clens( C )
  N = length(C)
  K = N
  for i = 3:N
    for k = 1:min(K,i)
      @printf("i:%3d  k:%3d  cld(i,k):%3d  i-k+1:%3d  diff:%3d  len:%3d\n",i,k,cld(i,k),i-k+1,(i-k+1)-cld(i,k)+1,length(C[i][k]))
    end
  end
end

function iters(N::Int64, K::Int64 )
  for i = 3:N
    for k = 1:min(K,i)
      #println("diff: ",(i-k+1) - cld(i,k) + 1 )
      @printf("i:%3d  k:%3d  cld(i,k):%3d  i-k+1:%3d  diff:%3d\n",i,k,cld(i,k),i-k+1,(i-k+1) - cld(i,k) + 1)
      #@printf("i:%3d  k:%3d  cld(i,k):%3d  i-k+1:%3d\n",i,k,cld(i,k),i-k+1)
      #for j = cld(i,k):(i-k+1)
        #@printf("%3d %3d %3d\n",i,k,j)
      #end
    end
  end
end

function unitVecs( N::MyInt )
  e1 = vcat([1],zeros(MyInt,N-1))
  result =  Array{MyInt,1}[e1]
  for i = 2:N
    ei = deepcopy(result[i-1])
    ei[i-1] = 0
    ei[i] = 1
    push!(result,ei)
  end
  return result
end
  
=#
