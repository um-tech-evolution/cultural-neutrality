#=
Estimates of theta = N_e mu based on Watterson "The homozygosity test of neutrality" (1978) 
  and Zeng, Fu, Shi, Wu "Statistical Tests for Detecting Positive Selection . . . " (2006).

allele_freqs is the vector of allele frequecies of the infinite alleles model population.
=#

export watterson

function watterson( allele_freqs::Vector{Int64} )
  n = sum(allele_freqs)
  sum_a_sq = 0
  for a in allele_freqs
    sum_a_sq += a^2
  end
  Float64(sum_a_sq)/Float64(n^2)
end

function watterson2( allele_freqs::Vector{Int64} )
  n = sum(allele_freqs)
  a_n = 0
  for i = 1:(n-1)
    a_n += 1.0/i
  end
end
