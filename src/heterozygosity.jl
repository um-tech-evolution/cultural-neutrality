@doc """ function estimated_heterozygocity( N::Int64, mu::Float64 )

Returns the theoreticali expected heterozygocity for given values of N and mu assuming a haploid model.
"""
function estimated_heterozygosity( N::Int64, mu::Float64 )
  1.0/(1.0 + 2*mu*N) 
end

@doc """ function heterozygosity( p::Array{Int64,2}, gen::Int64 )

Compute the empirical "heterozygosity" of a generation (i. e., row) of the allele matrix p which is returned by inf_alleles().
A pair of haploid genotypes is heterozygous if the aleles differ, and homozygous if the alleles are the same
Return the average heterozygoosity over all pairs of individuals.
"""
function heterozygosity( p::Array{Int64,2}, gen::Int64 ) 
  N = size(p,2)
  h = 0
  for i = 1:N
    for j = 1:N
      if p[gen,i] == p[gen,j]
        h += 1
      end
    end
  end
  return float(h)/N^2
end

@doc """ function heterozygosities( p::Array{Int64,2} )
Return a vector of the heterozygosities of all generations of the allele matrix p.
"""
function heterozygosities( p::Array{Int64,2} )
  nrows = size(p,1)
  [heterozygosity(p,j) for j = 1:nrows]
end
