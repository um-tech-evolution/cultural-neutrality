#=
Simulate the infinite alleles model (which is the Wright-Fisher model with infinite alleles mutation).
This is a single locus model.  Haploidy is assumed---which means that genotypes are not in diploid pairs.
=#

using DataStructures
using DataFrames
#include("proportional.jl")

@doc """ function inf_alleles( N::Int64, mu::Float64, ngens::Int64, burn_in::Int64=200 )

Run the infinite alleles model as a simulation. This is a 1-locus model.
Returns an array of "populations" where each population is a vector of Int64s.  
(Population size can vary from generation to generation.)
The initial population consists of the list of the first N positive integers.
mu is the mutation rate.  Mutation creates a new integer not used before.
burn_in is a number of generations of "burn in": these populations are not part of the returned array.
ngens is the number of populations in the returned array.
TODO:  add descriptions of K and C arguments
"""
function inf_alleles( N::Int64, mu::Float64, ngens::Int64; copy_funct::Function=neutral_copy, burn_in::Int64=200,
    K::Int64=0, C::Float64=0.0 )
  result = Array[ collect(1:N) ]
  new_id = N+1
  for g = 2:(ngens+burn_in)
    if K == 0 && C == 0.0 
      push!( result, copy_funct( result[g-1], mu) )
    elseif K > 0 && C > 0.0 
      push!( result, copy_funct( result[g-1], mu, K, C ) )
    else
      error("illegal copy_funct")
    end
  end
  result[burn_in+1:end]
end

function neutral_copy( lst::Vector{Int64}, mu::Float64 )
  new_id = maximum(lst)+1   # TODO:  add to argument list since this method can lead to subtle bugs.
  N = length(lst)
  p = zeros(Int64,N)
  for i = 1:N
    if rand() < mu
      p[i] = new_id
      new_id += 1
    else
      p[i] = lst[rand(1:N)]
    end
  end
  return p
end

@doc """ function topKlist( lst::Vector{Int64}, K::Int64 )

Return the list of the K most frequent elements of lst.
""" 
function topKlist( lst::Vector{Int64}, K::Int64 )
  c = counter(Int64)
  for x in lst
    push!(c,x)
  end
  result = sort( unique(lst), by=x->c[x], rev=true )
  result[1:min(K,length(result))]
end

@doc """ function conformist_copy( lst::Vector{Int64}, mu::Float64, K::Int64, C::Float64 )

Create a new "population" (list of integers) using the conformist copying method described in
"Biases in cultural transmission shape the turnover of popular traits" by Acerbi and Bentley
in Evolution and Human Behavior 35 (2014) 228–236.
mu  is the mutation rate, i. e., probablity of a new integer
C   is the probability of a conformist copy from the top K list
1-C  is the probability of a ranomd copy
returns  a list of N integers.
"""
function conformist_copy( lst::Vector{Int64}, mu::Float64, K::Int64, C::Float64 )
  N = length(lst)
  new_id = maximum(lst)+1   # TODO:  add to argument list since this method can lead to subtle bugs.
  result = zeros(Int64,N)
  topK = topKlist( lst, K )
  for i = 1:N
    if rand() < mu  # mutate
      result[i] = new_id
      new_id += 1
    elseif rand() < C  # conformist copy
      if lst[i] in topK
        result[i] = lst[i]  # copy if in topK
      else
        result[i] = lst[rand(1:length(lst))]
      end
    else     # random copy
      result[i] = lst[rand(1:N)]
    end
  end
  result
end

@doc """ function anti_conformist_copy( lst::Vector{Int64}, mu::Float64, K::Int64, C::Float64 )

Create a new "population" (list of integers) using the anti-conformist copying method described in
"Biases in cultural transmission shape the turnover of popular traits" by Acerbi and Bentley
in Evolution and Human Behavior 35 (2014) 228–236.
mu  is the mutation rate, i. e., probablity of a new integer
C   is the probability of a conformist copy from the top K list
1-C  is the probability of a ranomd copy
returns  a list of N integers.
"""
function anti_conformist_copy( lst::Vector{Int64}, mu::Float64, K::Int64, C::Float64 )
  N = length(lst)
  new_id = maximum(lst)+1   # TODO:  add to argument list since this method can lead to subtle bugs.
  result = zeros(Int64,N)
  topK = topKlist( lst, K )
  for i = 1:N
    if rand() < mu  # mutate
      result[i] = new_id
      new_id += 1
    elseif rand() < C  # anti_conformist copy
      if lst[i] in topK
        result[i] = lst[rand(1:length(lst))]
      else
        result[i] = lst[i]   # copy if not in topK
      end
    else     # random copy
      result[i] = lst[rand(1:N)]
    end
  end
  result
end


# not yet debugged or tested
function inf_alleles_with_selection( N::Int64, mu::Float64, ngens::Int64, 
      fit_funct::Function=neutral_fitness,  burn_in::Int64=200 )
  p = zeros(Int64, ngens+burn_in, N )
  p[1,:] = collect(1:N)
  new_id = N+1
  for g = 2:(ngens+burn_in)
    for i = 1:N
      if rand() < mu
        p[g,i] = new_id
        new_id += 1
      else
        p[g,i] = p[g-1,rand(1:N)]
      end
    end
  end
  p[burn_in+1:end,:]
end

@doc """ function neutral_fitness( x::Int64 )
  Every integer has fitness 1.0
"""
function neutral_fitness( x::Int64 )
  return 1.0
end

