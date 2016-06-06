#=
Estimates of theta = N_e mu based on Watterson "The homozygosity test of neutrality" (1978).

allele_freqs is the vector of allele frequecies of the infinite alleles model population.
=#

export watterson_homozygosity, watterson_theta

function watterson_homozygosity( allele_freqs::Vector{Int64} )
  n = sum(allele_freqs)
  sum_a_sq = 0
  for a in allele_freqs
    sum_a_sq += a^2
  end
  Float64(sum_a_sq)/Float64(n^2)
end

function watterson_theta( allele_freqs::Vector{Int64} )
  1.0/watterson_homozygosity( allele_freqs ) - 1.0
end