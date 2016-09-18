# Implements Frank M. Stewart's algorithm to generate allele count samples from the infinite alleles model
#    of a given size N and a given number of alleles K.  This is based on Warren Ewens' result
#    that the allele count configuration is independent of theta = 2 N mu (for haploid) or 4 N mu
#    (for diploid) where mu is the per locus mutation rate.
# A crucial function in the Stewart algorithm is B(K,N).  Julia functions whose name start with "B"
#    compute this function or a table of function values based on a formula or the recurrence 
#    relation given in his algorithm description.  "BS" functions return a single value, while
#    "BT" fucntions return a table of values.  
# Appendix to "STATISTICAL STUDIES ON PROTEIN POLYMORPHISM IN NATURAL POPULATIONS I. DISTRIBUTION 
#     OF SINGLE LOCUS HETEROZYGOSITY" by PAUL A. FUERST, RANAJIT CHAKRABORTY AND MASATOSHI NE1
#     Genetics, 1977.

export BSs, BSrec, BT, BS, BTr, BSr, Pcfg, PCcfg, Rsample, number_orders, number_ordersF
  
@doc """ function BSs( K::Int64, N::Int64 )
Based on Stewart's equation between A1 and A2.
Returns a rational number.
"""
function BSs( K::Int64, N::Int64 )
  if K > N 
    return 0 
  end
  return (-1)^(N-K)*factorial(K)*stirling1(N,K)//factorial(N)
end

@doc """ function BSr( K::Int64, N::Int64 )
Recursive version based on the recurrence given at the end of the appendix.
Returns a Rational
"""
function BSrec( K::Int64, N::Int64 )
  if K == 0 && N == 0
    return 1
  elseif K <= 0 || N <= 0 || K > N
    return 0
  end
  return (K*BSrec(K-1,N-1) + (N-1)*BSrec(K,N-1))//N
end

@doc """ function BT( K::Int64, N::Int64 )
Bottom up version also based on the recurrence.
Returns a table of of Float64's
"""
function BT( K::Int64, N::Int64 )
  #println("function BT  K: ", K, "  N: ",N)
  if K == 0 && N == 0
    return 1.0
  elseif K <= 0 || N <= 0 || K > N
    return 0.0
  end
  result = zeros(Float64,K,N)
  result[1,1] = 1.0
  for k = 1:K
    result[k,k] = 1//1
    for n = k+1:N
      if k > 1 && n > 1
        result[k,n] = (k*result[k-1,n-1] + (n-1)*result[k,n-1])/n
      elseif n > 1
        result[k,n] = (n-1)*result[k,n-1]/n
      else
        result[k,n] = 0.0
      end
    end
  end
  result
end

@doc """ function BS( K::Int64, N::Int64 )
  Calls BT to build the table and then returns a Float64 from an entry in the table.
"""
function BS( K::Int64, N::Int64 )
  if K == 0 && N == 0
    return 1.0
  elseif K <= 0 || N <= 0 || K > N
    return 0.0
  end
  tbl = BT(K,N)
  return tbl == 0 ? 0 : tbl[K,N]
end

@doc """ function BTr( K::Int64, N::Int64 )
Bottom up version based on the recurrence.
Returns a table of Rational's 
"""
function BTr( K::Int64, N::Int64 )
  if K == 0 && N == 0
    return 1
  elseif K <= 0 || N <= 0 || K > N
    return 0
  end
  result = zeros(Rational{Int128},K,N)
  result[1,1] = 1//1
  for k = 1:K
    result[k,k] = 1//1
    for n = k+1:N
      if k > 1 && n > 1
        result[k,n] = (k*result[k-1,n-1] + (n-1)*result[k,n-1])//n
      elseif n > 1
        result[k,n] = (n-1)*result[k,n-1]//n
      else
        result[k,n] = 0//1
      end
    end
  end
  result
end

@doc """ function BSr( K::Int64, N::Int64 )
  Calls BTr to build the table and then returns a Rational from an entry in the table.
"""
function BSr( K::Int64, N::Int64 )
  if K == 0 && N == 0
    return 1
  elseif K <= 0 || N <= 0 || K > N
    return 0
  end
  tbl = BTr(K,N)
  return tbl == 0 ? 0 : tbl[K,N]
end

@doc """ function Pcfg( N::Int64, K::Int64, allele_counts::Config )
Returns the probability of an ordered allele count configuration whose length
    may be less than K and whose sum may be less than N.
Ordered means that allele_counts is not assumed to be in decreasing order of counts.
The function sprob from ewens.jl returns the probability of an unordered allele_counts.
Thus:  
   number_orders(allele_counts)*Pcfg( sum(allele_counts), length(allele_counts), allele_counts ) 
    == Pcfg( sum(allele_counts), length(allele_counts), allele_counts )
Returns a rational.
Inefficient because it computes two tables.  Other versions are more efficient.
Example:
julia> Pcfg( 6, 3, [2,3,1] )
4//45
"""
function Pcfg( N::Int64, K::Int64, allele_counts::Config )
  L = length( allele_counts )
  S = sum( allele_counts )
  P = prod( allele_counts )
  return BSr(K-L,N-S)//P//BSr(K,N)
end

@doc """ function Pcfg( N::Int64, K::Int64, allele_counts::Config, Btbl::Array{Float64,2} )
Returns the probability of an ordered allele count configuration whose length
    may be less than K and whose sum may be less than N.
Ordered means that different orders of allele_counts are counted separately.
Btbl would normally be computed by the function BT.
"""
function Pcfg( N::Int64, K::Int64, allele_counts::Config, Btbl::Array{Float64,2} )
  L = length( allele_counts )
  S = sum( allele_counts )
  P = prod( allele_counts )
  if K-L <= 0 || N-S <= 0
    num = K-L==0 && N-S==0 ? 1. : 0.0
  else
    num = Btbl[K-L,N-S]
  end
  return num/P/Btbl[K,N]
end

@doc """ function Pcfg( N::Int64, K::Int64, allele_counts::Vector{Int64}, Btbl::Array{Float64,2} )
Returns the probability of an ordered allele count configuration whose length
    may be less than K and whose sum may be less than N.
Ordered means that different orders of allele_counts are counted separately.
Btbl would normally be computed by the function BT.
"""
function Pcfg( N::Int64, K::Int64, allele_counts::Vector{Int64}, Btbl::Array{Float64,2} )
  L = length( allele_counts )
  S = sum( allele_counts )
  P = prod( allele_counts )
  if K-L <= 0 || N-S <= 0
    num = K-L==0 && N-S==0 ? 1. : 0.0
  else
    num = Btbl[K-L,N-S]
  end
  return num/P/Btbl[K,N]
end

@doc """ function Pcfg( allele_counts::Config, Btbl::Array{Float64,2} )
Assumes that N and K are the sum and length of allele_counts respectively.
"""
function Pcfg( allele_counts::Config, Btbl::Array{Float64,2} )
  N = convert(Int64,sum(allele_counts))
  K = length(allele_counts)
  Pcfg(N,K,allele_counts,Btbl)
end

@doc """ function Pcfg( allele_counts::Vector{Int64}, Btbl::Array{Float64,2} )
Assumes that N and K are the sum and length of allele_counts respectively.
"""
function Pcfg( allele_counts::Vector{Int64}, Btbl::Array{Float64,2} )
  N = convert(Int64,sum(allele_counts))
  K = length(allele_counts)
  Pcfg(N,K,allele_counts,Btbl)
end

@doc """ function Pcfg( N::Int64, K::Int64, allele_counts::Config, Btbl::Array{Rational{Int128},2} )
Returns the probability of an ordered allele count configuration whose length
    may be less than K and whose sum may be less than N.
Ordered means that different orders of allele_counts are counted separately.
Btbl would normally be computed by the function BTr.
"""
function Pcfg( N::Int64, K::Int64, allele_counts::Config, Btbl::Array{Rational{Int128},2} )
  L = length( allele_counts )
  S = sum( allele_counts )
  P = prod( allele_counts )
  if K-L <= 0 || N-S <= 0
    num = K-L==0 && N-S==0 ? 1 : 0
  else
    num = Btbl[K-L,N-S]
  end
  return num//P//Btbl[K,N]
end

@doc """ function PCcfg( N::Int64, K::Int64, allele_counts::Config )
Returns the conditional probability of allele_counts[end] given allele_counts[1:(end-1)].
Note that allele_counts is ordered which means that allele_counts is not assumed to be 
   in decreasing order of counts.
Inefficient because it computes 2 tables.
Returns a rational.
"""
function PCcfg( N::Int64, K::Int64, allele_counts::Config )
  L = length( allele_counts )
  S = sum( allele_counts )
  Sr = S - allele_counts[end]  # sum except for last
  tmp = BSr(K-L+1,N-Sr)
  return tmp > 0 ? BSr(K-L,N-S)//allele_counts[end]//tmp : 0
end

@doc """ function PCcfg( N::Int64, K::Int64, allele_counts::Config, Btbl::Array{Float64,2} )
Returns the conditional probability of allele_counts[end] given allele_counts[1:(end-1)]
Note that allele_counts is ordered which means that allele_counts is not assumed to be
   in decreasing order of counts.
The table Btbl would normally be computed by the function BT().
Returns a Float64.
"""
function PCcfg( N::Int64, K::Int64, allele_counts::Config, Btbl::Array{Float64,2} )
  #println("PCcfg: N: ",N,"  K: ",K,"  allele_counts: ",transpose(allele_counts))
  L = length( allele_counts )
  S = sum( allele_counts )
  Sr = S - allele_counts[end]  # sum except for last
  if K-L <= 0 || N-S <= 0
    return K-L==0 && N-S==0 ? 1.0 : 0.0
  end
  tmp = Btbl[K-L+1,N-Sr] 
  return tmp > 0.0 ? Btbl[K-L,N-S]/allele_counts[end]/tmp : 0.0
end

@doc """ function PCcfg( N::Int64, K::Int64, allele_counts::Config, Btbl::Array{Float64,2} )
Returns the conditional probability of allele_counts[end] given allele_counts[1:(end-1)]
Note that allele_counts is ordered which means that allele_counts is not assumed to be
   in decreasing order of counts.
The table Btbl would normally be computed by the function BT().
Returns a Float64.
"""
function PCcfg( N::Int64, K::Int64, allele_counts::Vector{Int64}, Btbl::Array{Float64,2} )
  #println("PCcfg: N: ",N,"  K: ",K,"  allele_counts: ",transpose(allele_counts))
  L = length( allele_counts )
  S = sum( allele_counts )
  Sr = S - allele_counts[end]  # sum except for last
  if K-L <= 0 || N-S <= 0
    return K-L==0 && N-S==0 ? 1.0 : 0.0
  end
  tmp = Btbl[K-L+1,N-Sr] 
  return tmp > 0.0 ? Btbl[K-L,N-S]/allele_counts[end]/tmp : 0.0
end

@doc """ PCcfg( N::Int64, K::Int64, allele_counts::Config, Btbl::Array{Rational{Int128},2} )
Returns the conditional probability of allele_counts[end] given allele_counts[1:(end-1)]
Note that allele_counts is ordered which means that allele_counts is not assumed to be
   in decreasing order of counts.
The table Btbl would normally be computed by the function BTr().
Returns a rational.
"""
function PCcfg( N::Int64, K::Int64, allele_counts::Config, Btbl::Array{Rational{Int128},2} )
  L = length( allele_counts )
  S = sum( allele_counts )
  Sr = S - allele_counts[end]  # sum except for last
  if K-L <= 0 || N-S <= 0
    return K-L==0 && N-S==0 ? 1 : 0
  end
  tmp = Btbl[K-L+1,N-Sr]
  return tmp > 0 ? Btbl[K-L,N-S]//allele_counts[end]//tmp : 0
end

@doc """ function function Rsample( N::Int64, K::Int64, Btbl::Array{Float64,2} )
This is top level of the implementation of Stewart's algorithm.
Returns a random ordered allele_counts configuration (not in decreasing order of counts).
The table Btbl would normally be computed by the function BT().
"""
function Rsample( N::Int64, K::Int64, Btbl::Array{Float64,2} )
  if K == 1
    return [N]
  end
  alleles = zeros(Int64,K)
  for k = 1:K
    r = rand()
    n = 1
    cumm = 0.0
    while n < N && r > cumm
      alleles[k] = n
      #println("k:",k,"  n:",n,"  alleles[1:k]:",alleles[1:k],"  type:",typeof(alleles[1:k]))
      pc = PCcfg( N, K, alleles[1:k], Btbl )
      cumm += pc
      #println("pc: ",pc,"  cumm: ",cumm)
      n += 1
    end
  end
  #println("Rsample alleles: ",transpose(alleles))
  alleles
end

@doc """ function factorial128( n::Int64 )
  Doesn't overflow up to n == 32
"""
function factorial128( n::Int64 )
  if n > 32
    error("overflow in factorial128")
  end
  f = Int128(1)
  for i = 2:n
    f *= i
  end
  f
end

@doc """ function factorialF( n::Int64 )
  Doesn't overflow up to n == 32
"""
function factorialF( n::Int64 )
  f = 1.0
  for i = 2:n
    f *= i
  end
  f
end

@doc """ function number_orders( allele_counts::Vector{ConfigInt} )
Returns the number of ordered allele_count configurations per unordered configurations.
Examples:
julia> number_orders([1,3,2])
6
julia> number_orders([2,2,2])
1
julia> number_orders([1,4,1])
3
TODO:  Revise algorithm to reduce level of overflow.
"""
function number_orders( allele_counts::Vector{ConfigInt} )
  counts = zeros(Int64,maximum(allele_counts))
  for a in allele_counts
    counts[a] += 1
  end
  #counts = count_duplicates( allele_counts )
  p = Int128(1)
  for c in counts
    if c > 1
      p *= factorial128(c)
    end
  end
  div(factorial128(length(allele_counts)),p)
end

@doc """ function number_ordersF( allele_counts::Vector{ConfigInt} )
Returns the number of ordered allele_count configurations per unordered configurations as a Float64
Examples:
julia> number_ordersF([1,3,2])
6
julia> number_ordersF([2,2,2])
1
julia> number_ordersF([1,4,1])
3
TODO:  Revise algorithm to reduce level of overflow.
"""
function number_ordersF( allele_counts::Vector{ConfigInt} )
  counts = zeros(Int64,maximum(allele_counts))
  for a in allele_counts
    counts[a] += 1
  end
  #counts = count_duplicates( allele_counts )
  p = 1.0
  for c in counts
    if c > 1
      p *= factorialF(c)
    end
  end
  factorialF(length(allele_counts))/p
end

# Move to test functions or delete
function test_Rsample( N::Int64, K::Int64, nreps::Int64 )
  Btbl = BT(K,N)
  dictc = Dict{Config,Int64}()
  dictf = Dict{Config,Float64}()
  for K = 1:nreps
    sample = sort(Rsample(N,K,Btbl),rev=true)
    val = get(dictc,sample,0)
    dictc[sample] = val+1
  end
  nrepsf = convert(Float64,nreps)
  ssum = 0.0
  for key in keys(dictc)
    tmp = dictc[key]/nrepsf
    dictf[key] = tmp
    ssum += tmp
  end
  #println("sum: ",ssum)   # should be close to 1.0
  dictp
end

# Move to test functions or delete
function test_montecarlo( allele_counts::Config, nreps::Int64 )
  #N = sum(allele_counts)
  dict = test_Rsample( sum(allele_counts), length(allele_counts), nreps )
  ac32 = Int32[ c for c in allele_counts ]
  sr = ewens_montecarlo( Int32(nreps), ac32 )
  pkey = dict[allele_counts]/number_orders(allele_counts)
  println("pkey: ",pkey)
  ssum = 0.0
  for k in keys(dict)
    if dict[k]/number_ordered(k) <= pkey
      ssum += dict[k]
    end
  end
  (allele_counts,pkey,ssum,sr.probability)
end
  
# Move to test functions or delete
@doc """ function Bsum3( N::Int64, K::Int64 )
Test function used in implementing Rsample.
Computes the cummulative sum of conditional configurations 3 levels deep.
All three cummulative sums should be approximately 1.0.
"""
function Bsum3( N::Int64, K::Int64 )
  Btbl = BT(K,N)
  if K == 1
    return [N]
  end
  alleles = zeros(Int64,3)
  cumm1 = 0
  cumm2 = 0
  cumm3 = 0
  for n1 = 1:N
    alleles[1] = n1
    pc1 = PCcfg( N, K, alleles[1:1], Btbl )
    if K > 1 && pc1 > 0  
      cumm1 += pc1
      for n2 = 1:N
        alleles[2] = n2
        pc2 = PCcfg( N, K, alleles[1:2], Btbl )
        #println("n1:",n1,"  n2:",n2,"  pc1:",pc1,"  pc2:",pc2)
        if K > 2 && pc2 > 0
          cumm2 += pc1*pc2
          for n3 = 1:N
            alleles[3] = n3
            pc3 = PCcfg( N, K, alleles[1:3], Btbl )
            #println("n1:",n1,"  n2:",n2,"  n3:",n3,"  pc1:",pc1,"  pc2:",pc2,"  pc3:",pc3)
            cumm3 += pc1*pc2*pc3
            if pc1*pc2*pc3 > 0
              println("n1:",n1,"  n2:",n2,"  n3:",n3,"  pc1*pc2*pc3:",pc1*pc2*pc3)
            end
          end
        end
      end
    end
  end
  println("cumm1: ",cumm1)
  println("cumm2: ",cumm2)
  println("cumm3: ",cumm3)
end
