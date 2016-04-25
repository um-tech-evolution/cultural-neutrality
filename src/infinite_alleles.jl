#=
Simulate the infinite alleles model (which is the Wright-Fisher model with infinite alleles mutation).
This is a single locus model.  Haploidy is assumed---which means that genotypes are not in diploid pairs.
=#

using DataStructures
using DataFrames

@doc """ function inf_alleles( N::Int64, mu::Float64, ngens::Int64, burn_in::Int64=200 )

Run the infinite alleles model as a simulation. This is a 1-locus model.
The result is a ngens by N array of allele values where rows are generations and columns are individuals.
N is the population size.  
mu is the mutation rate.
ngens is the number of generations.
burn_in is a number of generations of "burn in": these generations are not part of the returned matrix
   and are not counted in ngens.
"""
function inf_alleles( N::Int64, mu::Float64, ngens::Int64, burn_in::Int64=200 )
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

