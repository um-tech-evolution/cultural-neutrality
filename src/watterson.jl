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

# Note:  primary version in  ewens/ewens.jl
# The alpha_j of the Watterson and Slatkin papers
# TODO:  A dictionary might be more efficient than an array
function w_alpha( allele_freqs::Vector{Int64} )
  k = length(allele_freqs)
  n = sum(allele_freqs)
  alpha = zeros(Int64,n)
  i = 1
  while i <= k
    alpha[allele_freqs[i]] += 1
    j = 1
    while i+j <= k && allele_freqs[i+j] == allele_freqs[i]
      alpha[allele_freqs[i]] += 1
      j += 1
    end
    i = i+j
  end
  alpha
end
  
    
