# Simple simulation of repeated proportional selection
# Initial population has k individuals of higher selection coefficient s while
#    all other individuals have selection coefficient 1.
# Demonstrates that for strong selection (i. e, s=2) and for k=1, the probability
#    of fixation of the more fit individual almost does not depend on population size.

using DataStructures
include("../src/aliases.jl")
include("../src/freq_scaled_fitness.jl")

@doc """ s_poplist()
  N =popsize
  ngens = a limit on the number of generations.  Should be set to something very large.
  k = number of individuals with selection coefficient s (remaining have sel coeff 1)
  dfe = fitness function  
  alt = true means use propsel_alt as a check of correctness of propsel
Runs successive generations of proportional selection until fixation starting with an 
  initial population with k copies of 2 and N-k copies of 1.
  Returns a pair (u,gen) where u is the individual (1 or 2) that fixed, 
    and gen is the generation of fixation.
"""
function s_poplist( N::Int64, ngens::Int64, k::Int64, dfe::Function, alt::Bool=false )
  fitness_table = Dict{Int64,Float64}()
  new_pop = vcat([2 for i = 1:k],[ 1 for i = 1:(N-k)])
  poplist= Population[ new_pop ]
  for p in new_pop
    fit = dfe_fitness(p, dfe, fitness_table )  # set fitness of p to be dfe(p).
  end
  for g = 2:ngens
    if alt
      new_pop = propsel_alt( poplist[g-1], dfe )
    else
      new_pop = propsel( poplist[g-1], dfe, fitness_table )
    end
    Base.push!(poplist,new_pop)
    u = unique(new_pop)
    if length(u) == 1
      return (u[1],g-1)
    end
  end
  poplist
end

@doc """ function propsel_alt()
  Roulette wheel version of proportial selection
  Less efficient in almost all circumstances than the above version.
  Written as check of correctness for the above version
"""
# Also in freq_scaled_fitness.jl
function propsel_alt( pop::Population, dfe::Function)
  fitness_table = Dict{Int64,Float64}()
  N = length(pop)
  new_pop = zeros(Int64,N)
  fitdict = Dict{Int64,Float64}()
  popctr = pop_counter( pop )
  #println("popctr: ",popctr)
  for p in pop
    fit = dfe_fitness(p, dfe, fitness_table )  # set fitness of p to be dfe(p).
  end
  #println("fitness_table: ",fitness_table)
  sum_fitness = 0.0
  for p in keys(popctr)
    dfit = dfe_fitness( p, dfe, fitness_table )*popctr[p]
    #println("p: ",p,"  dfit: ",dfit)
    sum_fitness += dfit
  end
  #println("sum_fitness: ",sum_fitness)
  for p in keys(popctr)
    fitdict[p] = dfe_fitness( p, dfe, fitness_table )*popctr[p]/sum_fitness
  end
  #println("fitdict: ",fitdict)
  for i = 1:N
    r = rand()
    #println("i: ",i,"  r: ",r)
    sumfit = 0.0
    for p in keys(fitdict)
      sumfit += fitdict[p]
      #println("p: ",p," fitdict: ",fitdict[p],"  sumfit: ",sumfit)
      if r < sumfit
        new_pop[i] = p
        break
      end
    end
  end
  #println("new_pop: ",new_pop)
  new_pop
end

function pop_counter( pop::Population )
  c = Dict{Int64,Int64}()
  for x in pop
    c[x] = get( c, x, 0 ) + 1
  end
  c
end

@doc """ dfe_fitness()
  Copied from freq_scaled_fitness
"""
function dfe_fitness( p::Int64, dfe::Function, fitness_table::Dict{Int64,Float64} )
  val = get(fitness_table,p,-1.0)
  #println("val; ",val)
  if val == -1.0   # p is not in fitness table
    val = dfe( p )
    fitness_table[p] = val
  end
  #println("tbl: ",fitness_table)
  val
end

@doc """ run_trials()
  N =popsize
  k = number of individuals with selection coefficient s (remaining have sel coeff 1)
  T = number of trails
  s = sel coeff
  alt = true means use propsel_alt as a check of correctness of propsel
  fixed = true  means sel coeff of different individual is s
  fixed = false  means sel coeff of different individual is 1.0 + s/N

Example run:
julia> run_trials(128,1,10000,2.)
sel coefs: 1: 1.0 2: 2.0
fixes: 1: 0.2106  2: 0.7894
ngens: 1: 1.564102564102564  2: 14.901570813275907
"""
function run_trials( N::Int64, k::Int64, T::Int64, s::Float64=0.0, alt::Bool=false; 
    fixed::Bool=true )
  ngens = 100*N
  function dfe_s(x)
    if x == 2
      return  fixed ? s : 1.0+s/N
    else
      return 1.0
    end
  end
  println("sel coefs: 1: ",dfe_s(1)," 2: ",dfe_s(2))
  count_1 = 0
  count_2 = 0
  gens_1 = 0
  gens_2 = 0
  for t = 1:T
    pl = s_poplist(N,ngens,k,dfe_s,alt)
    count_1 += pl[1]==1?1:0
    count_2 += pl[1]==2?1:0
    gens_1 += pl[1]==1?pl[2]:0
    gens_2 += pl[1]==2?pl[2]:0
  end
  FT = Float64(T)
  println("fixes: 1: ",count_1/FT,"  2: ",count_2/FT)
  println("ngens: 1: ", gens_1/count_1,"  2: ",gens_2/count_2)
end

# Not really used
function dfe_simple( x:: Int64, N::Int64; fit_inc::Float64=1.0/N )
  if x == 1 
    return 1.0 + fit_inc 
  else
    return 1.0
  end
end

function dfex(x)
  return Float64(x)
end

