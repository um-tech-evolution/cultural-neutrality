using DataStructures
# TODO:  put these type aliases in a separate file

export topKlist, bottomKlist, topKset, bottomKset, turnover, power_conformist_poplist, acerbi_conformist_poplist, 
    acerbi_anti_conformist_poplist, acerbi_mixed_conformist_poplist, power_mixed_conformist_poplist


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

@doc """ function bottomKlist( pop::Population, K::Int64 )
Return the list of the K least frequent elements of pop.
""" 
function bottomKlist( pop::Population, K::Int64 )
  c = counter(Int64)
  for x in pop
    push!(c,x)
  end
  result = sort( unique(pop), by=x->c[x] )
  result[1:min(K,length(result))]
end

@doc """ function Kset( pop::Population, K::Int64 )
Return the setof the K most frequent elements of pop.
""" 
function topKset( pop::Population, K::Int64 )
  Set(topKlist( pop, K ))
end

@doc """ function bottomKset( pop::Population, K::Int64 )
Return the setof the K most frequent elements of pop.
""" 
function bottomKset( pop::Population, K::Int64 )
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
N_mu:   per-population mutation rate.  The per-individual mutation rate mu = N_mu/N.
ngens:  number of Populations in returned list
Return a Population list using the power conformist copying method.
If  conformist_power == 0.0,  neutral, i. e., neither conformity nor anti-conformity
    conformist_power > 0.0,  conformist:  more likely to choose more frequent alleles from the previous generation
    conformist_power < 0.0,  anti-conformist:  more likely to choose less frequent alleles from the previous generation
See the documentation and code of function freq_scaled_fitness() for more details.
If combine==true, returns a single combined population of the populations from all ngens generations after burn in.
If combine==false, returns a list of ngens populations.
"""

function power_conformist_poplist( N::Int64, N_mu::Float64, ngens::Int64; burn_in::Float64=1.0, conformist_power::Float64=0.0,
    uniform_start::Bool=false, combine::Bool=true )
  #println("conformist power: ",conformist_power)
  int_burn_in = Int(round(N*burn_in))
  mu = N_mu/N
  if uniform_start  # All allele values start with the same value.  Start with a selective sweep.
    poplist= Population[ Int64[1 for i = 1:N] ]
    new_id = 2
  else
    poplist= Population[ collect(1:N) ]
    new_id = N+1
  end
  if combine
    pop_result = Population()
  end
  w = fill(1.0,N)
  for g = 2:(ngens+int_burn_in)
    fitness = freq_scaled_fitness( poplist[g-1], w, conformist_power )
    new_pop = propsel( poplist[g-1], fitness )
    for i in 1:N
      if rand() < mu
        new_pop[i] = new_id
        new_id += 1
      end
    end
    push!( poplist, new_pop )
    if combine && g >= int_burn_in+1
      pop_result = vcat( pop_result, new_pop )
    end
  end
  if combine
    return pop_result
  else
    return poplist[int_burn_in+1:end]
  end
end

@doc """ function acerbi_conformist_poplist( N::Int64, mu::Float64, ngens::Int64, C::Float64; toplist_size::Int64, burn_in::Float64=1.0,
    uniform_start::Bool=false )

Return a Population list using the conformist copying method described in
"Biases in cultural transmission shape the turnover of popular traits" by Acerbi and Bentley
in Evolution and Human Behavior 35 (2014) 228–236.
mu  is the mutation rate, i. e., probablity of a new integer allele 
C   is the probability of a conformist copy from the toplist
1-C  is the probability of a random copy
toplist_size   is the size of the toplist to use (Acerbi uses 10 for this value).
"""

function acerbi_conformist_poplist( N::Int64, N_mu::Float64, ngens::Int64, acer_C::Float64; toplist_size::Int64=10, burn_in::Float64=1.0,
    uniform_start::Bool=false, combine::Bool=true )
  int_burn_in = Int(round(N*burn_in))
  mu = N_mu/N
  if uniform_start  # All allele values start with the same value.  Start with a selective sweep.
    poplist= Population[ Int64[1 for i = 1:N] ]
    new_id = 2
  else
    poplist= Population[ collect(1:N) ]
    new_id = N+1
  end
  if combine
    pop_result = Population()
  end
  for g = 2:(ngens+int_burn_in)
    result = zeros(Int64,N)
    toplist = topKlist( poplist[g-1], toplist_size )
    for i = 1:N
      if rand() < mu  # mutate
        result[i] = new_id
        new_id += 1
      elseif rand() < acer_C  # conformist copy
        if poplist[g-1][i] in toplist
          result[i] = poplist[g-1][i]  # copy if in toplist
        else
          result[i] = poplist[g-1][rand(1:N)]
        end
      else     # random copy
        result[i] = poplist[g-1][rand(1:N)]
      end
    end
    push!( poplist, result )
    if combine && g >= int_burn_in+1
      pop_result = vcat( pop_result, result )
    end
  end
  if combine
    return pop_result
  else
    return poplist[int_burn_in+1:end]
  end
end

@doc """ function acerbi_anti_conformist_poplist( N::Int64, N_mu::Float64, ngens::Int64, C::Float64; toplist_size::Int64=1, burn_in::Float64=1.0,
    uniform_start::Bool=false )
Return a Population list using the anti-conformist copying method described in
"Biases in cultural transmission shape the turnover of popular traits" by Acerbi and Bentley
in Evolution and Human Behavior 35 (2014) 228–236.
mu  is the mutation rate, i. e., probablity of a new integer allele 
C   is the probability of a random copy
1-C  is the probability of a non-copy:  they retain their trait
toplist_size   is the size of the toplist to use (Acerbi uses 10 for this value).
TODO: Implement what might be called bottomlist anticonformism----which is different from Acerbi anti-conformism.
"""

function acerbi_anti_conformist_poplist( N::Int64, N_mu::Float64, ngens::Int64, C::Float64; toplist_size::Int64=0, burn_in::Float64=1.0,
    uniform_start::Bool=false, combine::Bool=true )
  int_burn_in = Int(round(N*burn_in))
  mu = N_mu/N
  if uniform_start  # All allele values start with the same value.  Start with a selective sweep.
    poplist= Population[ Int64[1 for i = 1:N] ]
    new_id = 2
  else
    poplist= Population[ collect(1:N) ]
    new_id = N+1
  end
  if combine
    pop_result = Population()
  end
  for g = 2:(ngens+int_burn_in)
    result = zeros(Int64,N)
    toplist = topKlist( poplist[g-1], toplist_size )
    for i = 1:N
      if rand() < mu  # mutate
        result[i] = new_id
        new_id += 1
      elseif rand() < C  # anti_conformist copy
        if poplist[g-1][i] in toplist
          result[i] = poplist[g-1][rand(1:N)]
        else
          result[i] = poplist[g-1][i]   # retain trait from previous generation
        end
      else     # random copy
        result[i] = poplist[g-1][rand(1:N)]
      end
    end
    push!( poplist, result )
    if combine && g >= int_burn_in+1
      pop_result = vcat( pop_result, result )
    end
  end
  if combine
    return pop_result
  else
    return poplist[int_burn_in+1:end]
  end
end

@doc """ function bottomlist_anti_conformist_poplist( N::Int64, N_mu::Float64, ngens::Int64, C::Float64; bottomlist_size::Int64=1, burn_in::Float64=1.0,
    uniform_start::Bool=false )
Return a Population list using the anti-conformist copying method described in
"Biases in cultural transmission shape the turnover of popular traits" by Acerbi and Bentley
in Evolution and Human Behavior 35 (2014) 228–236.
mu  is the mutation rate, i. e., probablity of a new integer allele 
C   is the probability of a random copy
1-C  is the probability of a non-copy:  they retain their trait
toplist_size   is the size of the toplist to use (Acerbi uses 10 for this value).
TODO: Implement what might be called bottomlist anticonformism----which is different from Acerbi anti-conformism.
"""

function bottomlist_anti_conformist_poplist( N::Int64, N_mu::Float64, ngens::Int64, C::Float64; bottomlist_size::Int64=1, burn_in::Float64=1.0,
    uniform_start::Bool=false, combine::Bool=true )
  int_burn_in = Int(round(N*burn_in))
  mu = N_mu/N
  if uniform_start  # All allele values start with the same value.  Start with a selective sweep.
    poplist= Population[ Int64[1 for i = 1:N] ]
    new_id = 2
  else
    poplist= Population[ collect(1:N) ]
    new_id = N+1
  end
  if combine
    pop_result = Population()
  end
  for g = 2:(ngens+int_burn_in)
    result = zeros(Int64,N)
    bottomlist = bottomKlist( poplist[g-1], bottomlist_size )
    for i = 1:N
      if rand() < mu  # mutate
        result[i] = new_id
        new_id += 1
      elseif rand() < C  # anti_conformist copy
        if poplist[g-1][i] in bottomlist
          result[i] = poplist[g-1][i]   # retain trait from previous generation
        else
          result[i] = poplist[g-1][rand(1:N)]
        end
      else     # random copy
        result[i] = poplist[g-1][rand(1:N)]
      end
    end
    push!( poplist, result )
    if combine && g >= int_burn_in+1
      pop_result = vcat( pop_result, result )
    end
  end
  if combine
    return pop_result
  else
    return poplist[int_burn_in+1:end]
  end
end

# TODO:  Check if Kandler's model of conformity agrees with Acerbi's.
@doc """ function acerbi_mixed_conformist_poplist( N::Int64, N_mu::Float64, ngens::Int64, 
    conformist_prob::Float64, anti_conformist_prob::Float64; 
    toplist_size::Int64=1, bottomlist_size::Int64=1, bottom::Bool=true,
    burn_in::Float64=2.0, uniform_start::Bool=false, combine::Bool=true )
Return either a list of ngens populations, or a combined population depending on whether "combine" is true.
The population will be generated by a mix of random copying, Acerbi conformist copying, and 
bottomlist or Acerbi anti-conformist copying.
bottom==true means use a bottomlist for anti-conformism
bottom==false means use Acerbi's method of anti-conformism
n = sample size
N = popsize
N_mu = mutation rate per population, i. e.,  mu*N.  Thus, mu = N_mu/N.
"""
function acerbi_mixed_conformist_poplist( N::Int64, N_mu::Float64, ngens::Int64, 
    conformist_prob::Float64, anti_conformist_prob::Float64; 
    toplist_size::Int64=1, bottomlist_size::Int64=1, bottom::Bool=true,
    burn_in::Float64=2.0, uniform_start::Bool=false, combine::Bool=true )
  if conformist_prob + anti_conformist_prob > 1.0
    error("conformist_prob + anti_conformist_prob must be less than 1.0 in power_mixed_conformist.")
  end
  int_burn_in = Int(round(N*burn_in))
  mu = N_mu/N
  if uniform_start  # All allele values start with the same value.  Start with a selective sweep.
    poplist= Population[ Int64[1 for i = 1:N] ]
    new_id = 2
  else
    poplist= Population[ collect(1:N) ]
    new_id = N+1
  end
  if combine
    pop_result = Population()
  end
  for g = 2:(ngens+int_burn_in)
    toplist = topKlist( poplist[g-1], toplist_size )
    if bottom
      bottomlist = bottomKlist( poplist[g-1], bottomlist_size )
    end
    result = zeros(Int64,N)
    for i = 1:N
      rnd = rand()
      if rnd < mu  # mutate
        result[i] = new_id
        new_id += 1
      elseif rnd < mu+conformist_prob  # conformist copy
        if poplist[g-1][i] in toplist
          result[i] = poplist[g-1][i]  # copy if in toplist
        else
          result[i] = poplist[g-1][rand(1:N)]
        end
      elseif rnd < mu+conformist_prob+anti_conformist_prob
        if bottom
          if poplist[g-1][i] in bottomlist
            result[i] = poplist[g-1][i]   # retain trait from previous generation
          else
            result[i] = poplist[g-1][rand(1:N)]
          end
        else
          if poplist[g-1][i] in toplist
            result[i] = poplist[g-1][rand(1:length(poplist[g-1]))]
          else
            result[i] = poplist[g-1][i]   # retain trait from previous generation
          end
        end
      else     # random copy
        result[i] = poplist[g-1][rand(1:N)]
      end
    end
    push!( poplist, result )
    if combine && g >= int_burn_in+1
      pop_result = vcat( pop_result, result )
    end
  end
  poplist[int_burn_in+1:end]
  if combine
    return pop_result
  else
    return poplist[int_burn_in+1:end]
  end
end

@doc """ function power_mixed_conformist_poplist( N::Int64, N_mu::Float64, ngens::Int64, 
    conformist_prob::Float64, anti_conformist_prob::Float64; 
    conformist_power::Float64=0.0, anti_conformist_power::Float64=0.0,
    burn_in::Float64=2.0, uniform_start::Bool=false, combine::Bool=true )
Return either a list of ngens populations, or a combined population depending on whether "combine" is true.
The population will be generated by a mix of random copying, power conformist copying, and 
power anti-conformist copying.
"""
function power_mixed_conformist_poplist( N::Int64, N_mu::Float64, ngens::Int64, 
    conformist_prob::Float64, anti_conformist_prob::Float64; 
    conformist_power::Float64=0.0, anti_conformist_power::Float64=0.0,
    burn_in::Float64=2.0, uniform_start::Bool=false, combine::Bool=true )
  if conformist_power < 0.0
    error("conformist_power in power_mixed_conformist should be nonnegative.")
  end
  if anti_conformist_power > 0.0
    error("conformist_power in power_mixed_conformist should be nonpostive.")
  end
  if conformist_prob + anti_conformist_prob > 1.0
    error("conformist_prob + anti_conformist_prob must be less than 1.0 in power_mixed_conformist.")
  end
  int_burn_in = Int(round(N*burn_in))
  mu = N_mu/N
  if uniform_start  # All allele values start with the same value.  Start with a selective sweep.
    poplist= Population[ Int64[1 for i = 1:N] ]
    new_id = 2
  else
    poplist= Population[ collect(1:N) ]
    new_id = N+1
  end
  if combine
    pop_result = Population()
  end
  for g = 2:(ngens+int_burn_in)
    w_conf = fill(1.0,N)
    w_anti_conf = fill(1.0,N)
    rp = randperm(N)
    if conformist_prob > 0.0
      conf_fitness = freq_scaled_fitness( poplist[g-1], w_conf, conformist_power )
      new_conf_pop = propsel( poplist[g-1], conf_fitness )[rp]
    end
    if anti_conformist_prob > 0.0
      anti_conf_fitness = freq_scaled_fitness( poplist[g-1], w_anti_conf, anti_conformist_power )
      new_anti_conf_pop = propsel( poplist[g-1], anti_conf_fitness )[rp]
    end
    result = zeros(Int64,N)
    for i = 1:N
      rnd = rand()
      if rnd < mu  # mutate
        result[i] = new_id
        new_id += 1
      elseif rnd < mu+conformist_prob  # conformist copy (without replacement)
        result[i] = new_conf_pop[i]
      elseif rnd < mu+conformist_prob+anti_conformist_prob
        result[i] = new_anti_conf_pop[i]
      else     # random copy
        result[i] = poplist[g-1][rand(1:N)]
      end
    end
    push!( poplist, result )
    if combine && g >= int_burn_in+1
      pop_result = vcat( pop_result, result )
    end
  end
  poplist[int_burn_in+1:end]
  if combine
    return pop_result
  else
    return poplist[int_burn_in+1:end]
  end
end

