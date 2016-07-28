# Compute Stirling numbers of the first kind based on the recurrence relation given in
#   the Wikipedia article: https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind

export stirling1_rec, stirling1

@doc """ function stirling1_rec( n::Int64, k::Int64 )
Recursive algorithm to compute the Stirling numbers of the first kind based on the
   formula given in the above Wikipedia article.
Inefficient due to repeated calculations of the function with the same parameters.    
Also, overflows (without warning) on moderately large values of N.
This can be handled by using BigInts as in the non-recursive version.
"""
function stirling1_rec( n::Int64, k::Int64 )
  println("n: ",n,"  k: ",k)
  if k == 0 
    if n == 0
      return 1
    else
      return 0
    end
  end
  if n == 0
    return 0
  end
  return -(n-1)*stirling1_rec(n-1,k) + stirling1_rec(n-1,k-1)
end

#=  TODO:  delete these versions.
@doc """ function stirling1( N::Int64, K::Int64 )
Bottom up computation of the Stirling numbers of the first kind.
   This makes the above recursive algorithm efficient by avoiding repeated calculations.
The algorithm uses both O(N,K) space and time.
This function will integer overflow for fairly small values of N and K.
"""
function stirling1( N::Int64, K::Int64 )
  if K > N || N < 0 || K < 0
    return 0
  end
  if K == 0 
    if N == 0
      return 1
    else
      return 0
    end
  end
  if N == 0
    return 0
  end
  r1 = zeros(Int64,2)
  r1[1] = 1
  result = Array{Int64,1}[ r1 ]
  println("result: ",result)
  for n = 2:N
    r = zeros(Int64, min(K,n+1) )
    for k = 1:K
      println("n: ",n,"  k: ",k)
      t1 = k > 1 && k <= n ? result[n-1][k-1] : 0
      println("rtmp: ",rtmp)
      t2 = k <= n-1 ? result[n-1][k] : 0
      println("t2: ",t2)
      r[k] = -(n-1)*t2 + t1
    end
    println("r: ",r)
    push!(result,r)
  end
  result
end
=#  

#=
@doc """ function stirling1( N::Int64, K::Int64 )
Bottom up computation of the Stirling numbers of the first kind.
   This makes the above recursive algorithm efficient by avoiding repeated calculations.
The algorithm uses both O(N,K) space and time 
  (not counting the time needed to handle larger and larger BigInt's).
"""
function stirling1( N::Int64, K::Int64 )
  if N == 0 && K == 0
    return 1
  elseif K > N || N <= 0 || K <= 0
    return 0
  end
  result = zeros(Int128,N,K)
  result[1,1] = 1
  for n = 2:N
    for k = 1:K
      tmp = k > 1 ? result[n-1,k-1] : 0
      result[n,k] = -(n-1)*result[n-1,k] + tmp
    end
  end
  result
end
=#
  
@doc """ function stirling1( N::Int64, K::Int64 )
Bottom up computation of the Stirling numbers of the first kind.
   This makes the above recursive algorithm efficient by avoiding repeated calculations.
The algorithm uses both O(N,K) space and time 
  (not counting the time needed to handle larger and larger BigInts).
"""
function stirling1( N::Integer, K::Integer )
  if N == 0 && K == 0
    return 1
  elseif K > N || N <= 0 || K <= 0
    return 0
  end
  result = zeros(BigInt,N,K)
  result[1,1] = 1
  for n = 2:N
    for k = 1:K
      tmp = k > 1 ? result[n-1,k-1] : 0
      result[n,k] = -(n-1)*result[n-1,k] + tmp
    end
  end
  result[N,K]
end
  

