#using Combinatorics

# Example runs
#julia> sprob(c9)
#6.657599063709085454484015981523344731039365639770995967956985336753368741547127e-02
#
#julia> spl = sprob_list(ncfgs(8,3))
#5-element Array{BigFloat,1}:
# 8.528784648187631530812706036917000618643271207736826073713067316478830338105348e-02
# 9.594882729211087420042643923240938166311300639658848614072494669509594882729187e-02
# 2.558635394456289978678038379530916844349680170575692963752665245202558635394464e-01
# 3.070362473347547974413646055437100213219616204690831556503198294243070362473331e-01
# 2.558635394456289978678038379530916844349680170575692963752665245202558635394464e-01
#
#julia> sum(spl)
#9.999999999999999826855257810514727780414433730581784952787084983247030155345757e-01


export sprob, sprobr, ewens, ewensr
# 
#include("stirling1.jl")   # not needed when run with ../src/NeutralCulturalEvolution loaded

# Slatkin formula for the exact probability of an allele_count configuraiton
#   (based on Ewens) from Slatkin Exact paper.
# Note that these are probabilities given K which is the number of alleles (richness)
# Float version
function sprob( allele_counts::Vector{Int64} )
  K = length( allele_counts )
  n = sum(allele_counts)
  alpha = w_alpha( allele_counts )
  #println("alpha: ",alpha)
  result = 1.0
  for i = 1:n
    result *= (n-i+1)/i^alpha[i]/factorial(alpha[i])
  end
  result/abs(stirling1(n,K))
end

# Returns the list of sprob values for ConfigList cl which can be produced by ncfgs(N,K) 
#    where function ncfgs() is defined in slatpart.jl
# Should add to 1 if cl includes configs corresponding to specific values of N and of K
# Does not add to 1 if cl includes configs corresponding to multiple K
function sprob_list(cl::ConfigList)
  map(sprob,convert_config_list_to_Int64(cl))
end

# rational version

function sprobr( allele_counts::Vector{Int64} )
  k = length( allele_counts )
  n = sum(allele_counts)
  alpha = w_alpha( allele_counts )
  #println("alpha: ",alpha)
  result = 1//1
  for i = 1:n
    result *= (n-i+1)//i^alpha[i]//factorial(alpha[i])
  end
  result//abs(stirling1(n,k))
end

# Test confirgurations from Slatkin's Exact Correction paper
c9 = [4,4,3,2,1,1,1];
c27 = [9,2,1,1,1,1,1];
keith = [52,9,8,4,4,2,2,1,1,1,1,1,1,1,1];

# The alpha_j of the Watterson and Slatkin papers
# TODO:  A dictionary might be more efficient than an array
function w_alpha( allele_counts::Vector{Int64} )
  k = length(allele_counts)
  n = sum(allele_counts)
  alpha = zeros(Int64,n)
  i = 1
  while i <= k
    alpha[allele_counts[i]] += 1
    j = 1
    while i+j <= k && allele_counts[i+j] == allele_counts[i]
      alpha[allele_counts[i]] += 1
      j += 1
    end
    i = i+j
  end
  alpha
end


# Ewens sampling formula as given by https://en.wikipedia.org/wiki/Ewens%27s_sampling_formula
#   or  slide 118 of http://www.stats.ox.ac.uk/~didelot/popgen/lecture7.pdf

# Float version
function ewens( allele_counts, theta )
n = length(allele_counts)
prod1 = 1
for i = 0:(n-1)
  prod1 *= (theta+i)
end
prod2 = 1
for j = 1:n 
  prod2 *= theta^allele_counts[j]/factorial(allele_counts[j])/j^allele_counts[j]
end
factorial(n)*prod2/prod1
end

# Rational number version
function ewensr( allele_counts, theta )
n = length(allele_counts)
prod1 = 1
for i = 0:(n-1)
  prod1 *= (theta+i)
end
prod2 = 1
for j = 1:n 
  prod2 *= theta^allele_counts[j]//factorial(allele_counts[j])//j^allele_counts[j]
end
factorial(n)*prod2//prod1
end

A4 = ((4,0,0,0),(2,1,0,0),(0,2,0,0),(1,0,1,0),(0,0,0,1))
A5 = ((5,0,0,0,0),(3,1,0,0,0),(2,0,1,0,0),(1,2,0,0,0),(1,0,0,1,0),(0,1,1,0,0),(0,0,0,0,1))


#= Not used
function isum( lst )
  sum = 0
  for i = 1:length(lst)
    sum += i*lst[i]
  end
  sum
end
=#

#= Not used
function P2(alleles,theta)
n = length(alleles)
prod = 1
for i = 1:n
  p1 = factorial(alleles[i])//(i^alleles[i])
  println("i:",i,"  f:",factorial(alleles[i]),"  a:",i^alleles[i],"  p1:",p1)
  p2 = theta^alleles[i]
  println("p2:",p2)
  prod *= p1*p2
end
prod
end
=#

