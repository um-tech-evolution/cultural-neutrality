#=
Estimates of theta = N_e mu based on Watterson "The homozygosity test of neutrality" (1978).

allele_freqs is the vector of allele frequecies of the infinite alleles model population.
=#

export watterson_homozygosity, p_homozygosity, ns_homozygosity, watterson_theta, IQV,
    p_1_4, p_1_6, p_2_6, p_3_0, w_homoz, 
    p_homoz_1_4, p_homoz_1_6, p_homoz_2_6, p_homoz_3_0 


function watterson_homozygosity( allele_freqs::Config )
  watterson_homozygosity( map(Int64, allele_freqs ) )
end

function watterson_homozygosity( allele_freqs::Population )
  if length(allele_freqs) == 0
    return 0.0   # Not sure that this is valid
  end
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

@doc """ function ns_homozygosity()
Returns n_homozygosity computed over pop members with fitness 1.0, and
s_homozygosity computed over pop members with fitness != 1.0.
"""
function ns_homozygosity( pop::Population, dfe::Function, fitness_table::Dict{Int64,Float64})
  println("ns_homozygosity")
  n_pop = Population()
  s_pop = Population()
  for x in pop
    if dfe_fitness(x, dfe, fitness_table ) == 1.0
      Base.push!(n_pop,x)
    else
      Base.push!(s_pop,x)
    end
  end
  println("length(n_pop): ",length(n_pop),"  length(s_pop): ",length(s_pop))
  if length( n_pop ) > 0
    n_allele_freqs = pop_counts64( n_pop )
    n_n = sum(n_allele_freqs)
    sum_n_sq = 0.0
    for a in n_allele_freqs
      sum_n_sq += Float64(a)^2
    end
    n_homoz = Float64(sum_n_sq)/Float64(n_n^2)
  else
    n_homoz = 0.0
  end
  if length( s_pop ) > 0
    s_allele_freqs = pop_counts64( s_pop )
    n_s = sum(s_allele_freqs)
    sum_s_sq = 0.0
    for a in s_allele_freqs
      sum_s_sq += Float64(a)^2
    end
    s_homoz = Float64(sum_s_sq)/Float64(n_s^2)
  else
    s_homoz = 0.0
  end
  ( n_homoz, s_homoz )
end

@doc """ function heterozygosity_ratio()
H_T/H_0 ratio given as formula 5.25 page 250 of Hartl & Clark 4th ed.
Derived from Kimura's 1983 book.
See notes/2_17_17heterozygostity_ratio_simulation.txt
and notes/2_18_17heterozygosity_ratio.docx
"""
function heterozygosity_ratio(N,s,haploid=true)
  P = haploid ? 1.0 : 2.0  #  1.0 for haploid, 2.0 for diploid
  2*(2*P*N*s-1+exp(-2*P*N*s))/(2*P*N*s*(1-exp(-2*P*N*s)))
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

const w_homoz = watterson_homozygosity
  
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

