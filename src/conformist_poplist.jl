using DataStructures
# TODO:  put these type aliases in a separate file

export topKlist, topKset, turnover, power_conformist_poplist, acerbi_conformist_poplist, acerbi_anti_conformist_poplist


@doc """ function topKlist( pop::Population, K::Int64 )

Return the list of the K most frequent elements of pop.
""" 
function topKlist( pop::Population, K::Int64 )
  c = counter(Int64)
  for x in pop
    push!(c,x)
  end
  result = sort( unique(pop), by=x->c[x], rev=true )
  result[1:min(K,length(result))]
end

@doc """ function topKset( pop::Population, K::Int64 )

Return the setof the K most frequent elements of pop.
""" 
function topKset( pop::Population, K::Int64 )
  Set(topKlist( pop, K ))
end

@doc """ function turnover( pop1::Population, pop2::Population, K::Int64 )

The number of alleles entering the toplist plus the number of alleles leaving the toplist
   in the transition from pop1 to pop1.
This is the defintion of Evans and Giametto rather than the definition of Bentley.
Usually, this value is twice Bentley's value.
The exception is when the toplist has less than K elements and an allele leaves without being replaced.
"""

function turnover( pop1::Population, pop2::Population, K::Int64 )
  toplist1 = topKset( pop1, K )
  toplist2 = topKset( pop2, K )
  length(setdiff( toplist1, toplist2 )) + length(setdiff( toplist2, toplist1 ))
end


@doc """funcction power_conformist_poplist( N::Int64, mu::Float64, ngens::Int64; burn_in::Int64=0, conformist_power::Float64=0.0 )

N:      population size
mu:     mutation rate
ngens:  number of Populations in returned list
Return a Population list using the power conformist copying method.
If  conformist_power == 0.0,  neutral, i. e., neither conformity nor anti-conformity
    conformist_power > 0.0,  conformist:  more likely to choose more frequent alleles from the previous generation
    conformist_power < 0.0,  anti-conformist:  more likely to choose less frequent alleles from the previous generation
See the documentation and code of function freq_scaled_fitness() for more details.
"""

function power_conformist_poplist( N::Int64, mu::Float64, ngens::Int64; burn_in::Int64=0, conformist_power::Float64=0.0,
    uniform_start::Bool=false )
  if uniform_start  # All allele values start with the same value.  Start with a selective sweep.
    poplist= Population[ Int64[1 for i = 1:N] ]
    new_id = 2
  else
    poplist= Population[ collect(1:N) ]
    new_id = N+1
  end
  w = fill(1.0,N)
  for g = 2:(ngens+burn_in)
    fitness = freq_scaled_fitness( poplist[g-1], w, conformist_power )
    new_pop = propsel( poplist[g-1], fitness )
    for i in 1:N
      if rand() < mu
        new_pop[i] = new_id
        new_id += 1
      end
    end
    push!( poplist, new_pop )
  end
  poplist[burn_in+1:end]
end

@doc """ function acerbi_conformist_poplist( N::Int64, mu::Float64, ngens::Int64, K::Int64, C::Float64; burn_in::Int64=0 )

Return a Population list using the conformist copying method described in
"Biases in cultural transmission shape the turnover of popular traits" by Acerbi and Bentley
in Evolution and Human Behavior 35 (2014) 228–236.
mu  is the mutation rate, i. e., probablity of a new integer allele 
C   is the probability of a conformist copy from the top K list
1-C  is the probability of a random copy
K   is the size of the top K list to use.
"""

function acerbi_conformist_poplist( N::Int64, mu::Float64, ngens::Int64, K::Int64, C::Float64; burn_in::Int64=0,
    uniform_start::Bool=false )
  if uniform_start  # All allele values start with the same value.  Start with a selective sweep.
    poplist= Population[ Int64[1 for i = 1:N] ]
    new_id = 2
  else
    poplist= Population[ collect(1:N) ]
    new_id = N+1
  end
  for g = 2:(ngens+burn_in)
    result = zeros(Int64,N)
    topK = topKlist( poplist[g-1], K )
    for i = 1:N
      if rand() < mu  # mutate
        result[i] = new_id
        new_id += 1
      elseif rand() < C  # conformist copy
        if poplist[g-1][i] in topK
          result[i] = poplist[g-1][i]  # copy if in topK
        else
          result[i] = poplist[g-1][rand(1:N)]
        end
      else     # random copy
        result[i] = poplist[g-1][rand(1:N)]
      end
    end
    push!( poplist, result )
  end
  poplist[burn_in+1:end]
end

function acerbi_anti_conformist_poplist( N::Int64, mu::Float64, ngens::Int64, K::Int64, C::Float64; burn_in::Int64=0,
    uniform_start::Bool=false )
  if uniform_start  # All allele values start with the same value.  Start with a selective sweep.
    poplist= Population[ Int64[1 for i = 1:N] ]
    new_id = 2
  else
    poplist= Population[ collect(1:N) ]
    new_id = N+1
  end
  for g = 2:(ngens+burn_in)
    result = zeros(Int64,N)
    topK = topKlist( poplist[g-1], K )
    for i = 1:N
      if rand() < mu  # mutate
        result[i] = new_id
        new_id += 1
      elseif rand() < C  # anti_conformist copy
        if poplist[g-1][i] in topK
          result[i] = poplist[g-1][rand(1:length(poplist[g-1]))]
        else
          result[i] = poplist[g-1][i]   # copy if not in topK
        end
      else     # random copy
        result[i] = poplist[g-1][rand(1:N)]
      end
    end
    push!( poplist, result )
  end
  poplist[burn_in+1:end]
end

# TODO:  Check if Kandler's model of conformity agrees with Acerbi's.
