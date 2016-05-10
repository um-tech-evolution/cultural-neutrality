#=
Simulate the infinite alleles model (which is the Wright-Fisher model with infinite alleles mutation).
This is a single locus model.  Haploidy is assumed---which means that genotypes are not in diploid pairs.
=#

using DataStructures
using DataFrames
typealias Population Array{Int64,1}
typealias PopList Array{Array{Int64,1},1}

@doc """ function inf_alleles( N::Int64, mu::Float64, ngens::Int64, burn_in::Int64=200 )

Run the infinite alleles model as a simulation. This is a 1-locus model.
Allele values are represented by positive integers, and populations are list of positive integers.
The result is a list of Populations, one per generation.
N is the population size, or the initial population size if popsize_ratio != 1.0 (see below)  
mu is the mutation rate, or the initial mutation rate if popsize_ratio != 1.0 (see below)
ngens is the number of generations.
burn_in is a number of generations of "burn in": these generations are not part of the returned matrix
   and are not counted in ngens.
uniform_start == true means that the inital population is all ones.
popsize_ratio != 1.0 means that the popsize grows (or shrinks) geometrically with this ratio.
"""
function neutral_poplist( N::Int64, mu::Float64, ngens::Int64; burn_in::Int64=0, uniform_start::Bool=false,
    popsize_ratio::Float64=1.0 )
  if uniform_start  # All allele values start with the same value.  Start with a bottleneck.
    poplist= Population[ Int64[1 for i = 1:N] ]
    new_id = 2
  else
    poplist= Population[ collect(1:N) ]
    new_id = N+1
  end
  popsize = N
  for g = 2:(ngens+burn_in)
    previous_popsize = popsize
    popsize = (popsize_ratio == 1.0) ? N : max(1,Int(round(N*popsize_ratio^(g-1))))
    current_mu = (popsize_ratio == 1.0) ? mu : min(1.0, mu/popsize_ratio^(g-1))
    println("g: ",g,"  popsize: ",popsize,"  current mu: ",current_mu)
    result = zeros(Int64,popsize)
    for i = 1:popsize
      if rand() < current_mu
        result[i] = new_id
        new_id += 1
      else
        result[i] = poplist[g-1][rand(1:previous_popsize)]
      end
    end
    push!(poplist,result)
  end
  poplist[burn_in+1:end]
end


@doc """ function topKlist( pop::Population, K::Int64 )

Return the list of the K most frequent elements of pop.
""" 
function topKlist( pop::Population, K::Int64 )
  c = counter(Int64)
  for x in pop
    push!(c,x)
  end
  result = sort( unique(pop), by=x->c[x], rev=true )
  result[1:min(K,length(result))]
end

@doc """ function topKset( pop::Population, K::Int64 )

Return the setof the K most frequent elements of pop.
""" 
function topKset( pop::Population, K::Int64 )
  Set(topKlist( pop, K ))
end

@doc """ function turnover( pop1::Population, pop2::Population, K::Int64 )

The number of alleles entering the toplist plus the number of alleles leaving the toplist
   in the transition from pop1 to pop1.
This is the defintion of Evans and Giametto rather than the definition of Bentley.
Usually, this value is twice Bentley's value.
The exception is when the toplist has less than K elements and an allele leaves without being replaced.
"""

function turnover( pop1::Population, pop2::Population, K::Int64 )
  toplist1 = topKset( pop1, K )
  toplist2 = topKset( pop2, K )
  length(setdiff( toplist1, toplist2 )) + length(setdiff( toplist2, toplist1 ))
end

@doc """ function pop_counts( pop::Population )

Returns the sorted frequencies of the alleles of Population pop.
Example:  If pop = [5, 7, 9, 5, 4, 5, 7], then the returned list is [3, 2, 1, 1]
   because there are 3 5's, 2 7's, 1 9, and 1 4.  So the sum of the returned 
   list is the length of the population.
"""
function pop_counts( pop::Population )
  c = counter(Int64)
  for x in pop
    push!(c,x)
  end
  map( x->c[x], sort( unique(pop), by=x->c[x], rev=true ) )
end

@doc """ function poplist_counts( poplst::PopList )

Returns the sorted frequencies of the combined populations of poplist.
"""
function poplist_counts( poplst::PopList )
  combined_pop = Int64[]
  for pop in poplst
    combined_pop = vcat(combined_pop,pop)
  end
  pop_counts( combined_pop )
end
