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
N is the population size.  
mu is the mutation rate.
ngens is the number of generations.
burn_in is a number of generations of "burn in": these generations are not part of the returned matrix
   and are not counted in ngens.
"""
function neutral_poplist( N::Int64, mu::Float64, ngens::Int64; burn_in::Int64=0 )
  poplist= Population[ collect(1:N) ]
  new_id = N+1
  for g = 2:(ngens+burn_in)
    result = zeros(Int64,N)
    for i = 1:N
      if rand() < mu
        result[i] = new_id
        new_id += 1
      else
        result[i] = poplist[g-1][rand(1:N)]
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

function pop_counts( pop::Population )
  c = counter(Int64)
  for x in pop
    push!(c,x)
  end
  map( x->c[x], sort( unique(pop), by=x->c[x], rev=true ) )
end

function poplist_counts( poplst::PopList )
  combined_pop = Int64[]
  for pop in poplst
    combined_pop = vcat(combined_pop,pop)
  end
  pop_counts( combined_pop )
end
