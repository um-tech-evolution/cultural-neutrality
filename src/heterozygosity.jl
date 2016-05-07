@doc """ function estimated_heterozygocity( N::Int64, mu::Float64 )

Returns the theoreticali expected heterozygocity for given values of N and mu assuming a haploid model.
"""
function estimated_heterozygosity( N::Int64, mu::Float64 )
  1.0/(1.0 + 2*mu*N) 
end

@doc """ function heterozygosity( p::Array{Int64,2}, gen::Int64 )

Compute the empirical "heterozygosity" of a population of the array of populaitons returned by inf_alleles().
A pair of haploid genotypes is heterozygous if the aleles differ, and homozygous if the alleles are the same
Return the average heterozygoosity over all pairs of individuals.
"""
function heterozygosity( p, gen::Int64 ) 
  N = length(p[gen])
  h = 0
  for i = 1:N
    for j = 1:N
      if i != j && p[gen][i] == p[gen][j]
        h += 1
      end
    end
  end
  return float(h)/(N^2 - N)
end

@doc """ function heterozygosities( p::Array{Int64,2} )
Return a vector of the heterozygosities of all generations of the allele matrix p.
"""
function heterozygosities( p)
  nrows = length(p)
  [heterozygosity(p,j) for j = 1:nrows]
end
