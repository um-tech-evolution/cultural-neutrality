#=
Estimates of theta = N_e mu based on Watterson "The homozygosity test of neutrality" (1978).

allele_freqs is the vector of allele frequecies of the infinite alleles model population.
=#

export watterson_homozygosity, p_homozygosity, watterson_theta, IQV,
    p_1_4, p_1_6, p_2_6, p_3_0, w_homoz, 
    p_homoz_1_4, p_homoz_1_6, p_homoz_2_6, p_homoz_3_0 


function watterson_homozygosity( allele_freqs::Config )
  watterson_homozygosity( map(Int64, allele_freqs ) )
end

function watterson_homozygosity( allele_freqs::Population )
  n = sum(allele_freqs)
  sum_a_sq = 0
  for a in allele_freqs
    sum_a_sq += a^2
  end
  Float64(sum_a_sq)/Float64(n^2)
end

function watterson_theta( allele_freqs::Config )
  1.0/watterson_homozygosity( allele_freqs ) - 1.0
end

function watterson_theta( allele_freqs::Population )
  1.0/watterson_homozygosity( allele_freqs ) - 1.0
end

@doc """ function IQV()
Allan Wilcox's "index of quantitative variation" defined by formula (7) 
  of Mark E. Madsen's "Neutral Cultural Transmission in Time Averaged
Archaeological Assemblages" arXiv:1024.2043v2 (2012). 
A measure of "evenness".  IQV of a population with a single trait is 0,
while IQV of a population with an even distribution of traits is close to 1.0.
"""
function IQV( allele_freqs::Population )
  k = length(allele_freqs)
  if k == 1
    return 0.0
  end
  Fk = Float64(k)
  wh = watterson_homozygosity( allele_freqs )
  (1.0-wh)*Fk/(Fk-1.0)
end
  

function p_homozygosity( allele_freqs::Config, p::Float64 )
  p_homozygosity( map(Int64, allele_freqs ), p )
end

function p_homozygosity( allele_freqs::Population, p::Float64 )
  n = sum(allele_freqs)
  sum_a_sq = 0.0
  for a in allele_freqs
    sum_a_sq += Float64(a)^p
  end
  Float64(sum_a_sq)/Float64(n^p)
end

#=
function watterson_theta( allele_freqs::Config )
  1.0/watterson_homozygosity( allele_freqs ) - 1.0
end
=#

# Note:  primary version in  ewens/ewens.jl
# The alpha_j of the Watterson and Slatkin papers
# TODO:  A dictionary might be more efficient than an array
function w_alpha( allele_freqs::Config )
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

typealias w_homoz watterson_homozygosity
  
function p_1_4(cfg)
  p_homozygosity(cfg,1.4)
end

function p_1_6(cfg)
  p_homozygosity(cfg,1.6)
end

function p_2_6(cfg)
  p_homozygosity(cfg,2.6)
end

function p_3_0(cfg)
  p_homozygosity(cfg,3.0)
end

# for backward compatibility
function p_homoz_1_4(cfg)
  p_homozygosity(cfg,1.4)
end

function p_homoz_1_6(cfg)
  p_homozygosity(cfg,1.6)
end

function p_homoz_2_6(cfg)
  p_homozygosity(cfg,2.6)
end

function p_homoz_3_0(cfg)
  p_homozygosity(cfg,3.0)
end

