#=
Simulate the infinite alleles model (which is the Wright-Fisher model with infinite alleles mutation).
This is a single locus model.  Haploidy is assumed---which means that genotypes are not in diploid pairs.
=#
export neutral_poplist, pop_counts, poplist_counts

using DataStructures
#using DataFrames

@doc """ function neutral_poplist( N::Int64, mu::Float64, ngens::Int64; burn_in::Int64=N, uniform_start::Bool=false,
    popsize_ratio::Float64=1.0 )

Run the infinite alleles model as a simulation. This is a 1-locus model.
Allele values are represented by positive integers, and populations are list of positive integers.
The result is a list of Populations, one per generation.
N is the population size, or the initial population size if popsize_ratio != 1.0 (see below)  
mu is the mutation rate, or the initial mutation rate if popsize_ratio != 1.0 (see below)
ngens is the number of generations.
burn_in specifies the number of generations of "burn" in as a multiple of N (1.0 means 1.0*N generations of burn in): 
    these generations are not part of the returned matrix and are not counted in ngens.
uniform_start == true means that the inital population is all ones.
popsize_ratio != 1.0 means that the popsize grows (or shrinks) geometrically with this ratio.
"""
function neutral_poplist( N::Int64, mu::Float64, ngens::Int64; burn_in::Float64=1.0, uniform_start::Bool=false,
    popsize_ratio::Float64=1.0 )
  int_burn_in = Int(round(N*burn_in))
  if uniform_start  # All allele values start with the same value.  Start with a bottleneck.
    poplist= Population[ Int64[1 for i = 1:N] ]
    new_id = 2
  else
    poplist= Population[ collect(1:N) ]
    new_id = N+1
  end
  popsize = N
  for g = 2:(ngens+int_burn_in)
    previous_popsize = popsize
    popsize = (popsize_ratio == 1.0) ? N : max(1,Int(round(N*popsize_ratio^(g-1))))
    current_mu = (popsize_ratio == 1.0) ? mu : min(1.0, mu/popsize_ratio^(g-1))
    #println("g: ",g,"  popsize: ",popsize,"  current mu: ",current_mu)
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
  poplist[int_burn_in+1:end]
end

@doc """ function pop_counts( pop::Population )

Returns the sorted frequencies of the alleles of Population pop.
Example:  If pop = [5, 7, 9, 5, 4, 5, 7], then the returned list is [3, 2, 1, 1]
   because there are 3 5's, 2 7's, 1 9, and 1 4.  So the sum of the returned 
   list is the length of the population.
"""
function pop_counts( pop::Population )
  c = Dict{Int64,Int32}()
  for x in pop
    c[x] = get( c, x, 0 ) + 1
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
