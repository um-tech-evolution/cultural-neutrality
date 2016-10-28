#=
Simulate the infinite alleles model (which is the Wright-Fisher model with infinite alleles mutation).
This is a single locus model.  Haploidy is assumed---which means that genotypes are not in diploid pairs.
=#
export neutral_poplist, pop_counts8, pop_counts32,  pop_counts64, poplist_counts32, 
    poplist_counts64, simple_poplist, ewens_K_est, sample_population

using DataStructures
#using DataFrames

@doc """ function simple_poplist( N::Int64, mu_per_locus::Float64, ngens::Int64, burn_in::Int64 )
No-frills implmentation of Wright-Fisher infinite alleles model.
If combine==true, returns a single combined population of the populations from all ngens generations after burn in.
If combine==false, returns a list of ngens populations.
"""

function simple_poplist( N::Int64, mu_per_pop::Float64, ngens::Int64; burn_in::Float64=2.0, combine::Bool=true )
  int_burn_in = Int(round(N*burn_in))
  mu = mu_per_pop/N
  #println("int_burn_in: ",int_burn_in,"  mu: ",mu)
  poplist= Population[ collect(1:N) ]
  if combine
    pop_result = Population()
  end
  new_id = N+1
  for g = 2:(ngens+int_burn_in)
    result = zeros(Int64,N)
    for i = 1:N
      if rand() < mu
        result[i] = new_id
        new_id += 1
      else  # Choose a random element of the previous population
        result[i] = poplist[g-1][rand(1:N)]
      end
    end
    push!(poplist,result)
    if combine && g > int_burn_in+1
      println("combine: ", g)
      pop_result = vcat( pop_result, result )
    end
  end
  if combine
    return pop_result
  else
    return poplist[int_burn_in+1:end]
  end
end 

function simple_combined_pop( N::Int64, mu_per_pop::Float64, ngens::Int64; burn_in::Float64=1.0 )
  combine_pops( simple_poplist( N, mu_per_pop, ngens, burn_in=burn_in ) )
end

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

@doc """ function pop_counts32( pop::Population )
Returns the sorted frequencies of the alleles of Population pop.
Example:  If pop = [5, 7, 9, 5, 4, 5, 7], then the returned list is [3, 2, 1, 1]
   because there are 3 5's, 2 7's, 1 9, and 1 4.  So the sum of the returned 
   list is the length of the population.
"""
function pop_counts8( pop::Population )
  c = Dict{Int64,ConfigInt}()
  for x in pop
    c[x] = get( c, x, 0 ) + 1
  end
  map( x->c[x], sort( unique(pop), by=x->c[x], rev=true ) )
end

@doc """ function pop_counts32( pop::Population )
Returns the sorted frequencies of the alleles of Population pop.
Example:  If pop = [5, 7, 9, 5, 4, 5, 7], then the returned list is [3, 2, 1, 1]
   because there are 3 5's, 2 7's, 1 9, and 1 4.  So the sum of the returned 
   list is the length of the population.
"""
function pop_counts32( pop::Population )
  c = Dict{Int64,Int32}()
  for x in pop
    c[x] = get( c, x, 0 ) + 1
  end
  map( x->c[x], sort( unique(pop), by=x->c[x], rev=true ) )
end

@doc """ function pop_counts64( pop::Population )

Returns the sorted frequencies of the alleles of Population pop.
Example:  If pop = [5, 7, 9, 5, 4, 5, 7], then the returned list is [3, 2, 1, 1]
   because there are 3 5's, 2 7's, 1 9, and 1 4.  So the sum of the returned 
   list is the length of the population.
"""
function pop_counts64( pop::Population )
  c = Dict{Int64,Int64}()
  for x in pop
    c[x] = get( c, x, 0 ) + 1
  end
  map( x->c[x], sort( unique(pop), by=x->c[x], rev=true ) )
end

@doc """ function poplist_counts32( poplst::PopList )

Returns the sorted frequencies of the combined populations of poplist.
"""
function poplist_counts32( poplst::PopList )
  combined_pop = Int64[]
  for pop in poplst
    combined_pop = vcat(combined_pop,pop)
  end
  pop_counts32( combined_pop )
end

@doc """ function poplist_counts64( poplst::PopList )

Returns the sorted frequencies of the combined populations of poplist.
"""
function poplist_counts64( poplst::PopList )
  combined_pop = Int64[]
  for pop in poplst
    combined_pop = vcat(combined_pop,pop)
  end
  pop_counts64( combined_pop )
end

@doc """ function ewens_K_est( theta::Float64, N::Int64 )
Estimate K from theta and N using Ewens formula.
"""
function ewens_K_est( theta::Float64, N::Int64 )
  result = 1.0
  for i = 1:(N-1)
    result += theta/(theta+i)
  end
  result
end
  
@doc """ function sample_population( pop::Population, new_size::Int64 )
 Take a random sample (with replcement) of size new_size from a population
"""
function sample_population( pop::Population, new_size::Int64 )
  indices = rand(1:length(pop),new_size)
  new_pop = [ pop[i] for i in indices ]
end

@doc """ function combine_pops( poplist::Array{Vector{Int64},1} )
  Combine a list of populations into a single population
  Note:  Functions like poplist_counts64 may be preferable.
"""
function combine_pops( poplist::Array{Vector{Int64},1} )
  if length(poplist) == 1
    return poplist[1]
  end
  combined_pop = Int64[]
  for pop in poplst
    combined_pop = vcat(combined_pop,pop)
  end
  combined_pop
end
  
