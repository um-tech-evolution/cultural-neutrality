using DataStructures

export topKlist, bottomKlist, topKset, bottomKset, turnover, power_conformist_poplist, acerbi_conformist_poplist, 
    acerbi_anti_conformist_poplist, acerbi_mixed_conformist_poplist, power_mixed_conformist_poplist, 
    nearly_neutral_power_mixed_conformist_poplist


@doc """ function topKlist( pop::Population, K::Int64 )
Return the list of the K most frequent elements of pop.
If there are more than K of these, return the list of all of them.
""" 
function topKlist( pop::Population, K::Int64 )
  c = DataStructures.counter(Int64)
  for x in pop
    push!(c,x)
  end
  result = sort( pop, by=x->c[x], rev=true )
  j = K
  max = c[result[K]]
  while c[result[j]] == max && j < length(result)
    j += 1
  end
  result[1:j-1]
end

@doc """ function bottomKlist( pop::Population, K::Int64 )
Return the list of the K least frequent elements of pop.
If there are more than K of these, return the list of all of them.
""" 
function bottomKlist( pop::Population, K::Int64 )
  c = DataStructures.counter(Int64)
  for x in pop
    push!(c,x)
  end
  result = sort( pop, by=x->c[x] )
  j = K
  min = c[result[K]]
  while c[result[j]] == min && j < length(result)
    j += 1
  end
  result[1:j-1]
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

#=
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
=#

@doc """ function new_toplist_population( poplist::PopList, g::Int64, new_id::Vector{Int64},  N::Int64, mu::Float64,
      conformist_prob::Float64, anti_conformist_prob::Float64, toplist_size::Int64, bottomlist_size::Int64, bottom::Bool )
Generates a new population using the "choose from toplist" and "choose from bottomlist" methods 
    for conformity and anti-conformity.  These are used by Mesoudi-Lycett.
"""

function new_toplist_population( poplist::PopList, g::Int64, new_id::Vector{Int64},  N::Int64, mu::Float64,
      conformist_prob::Float64, anti_conformist_prob::Float64, toplist_size::Int64, bottomlist_size::Int64, bottom::Bool )
  toplist = topKlist( poplist[g-1], toplist_size )
  bottomlist = bottomKlist( poplist[g-1], bottomlist_size )
  result = zeros(Int64,N)
  for i = 1:N
    rnd = rand()
    if rnd < mu  # mutate
      result[i] = new_id[1]
      new_id[1] += 1
    elseif rnd < mu+conformist_prob  # conformist copy
      result[i] = toplist[rand(1:length(toplist))]
    elseif rnd < mu+conformist_prob+anti_conformist_prob  # anti-conformist copy
      result[i] = bottomlist[rand(1:length(toplist))]
    else  # random copy
      result[i] = poplist[g-1][rand(1:N)]
    end
  end
  return result
end

@doc """ function new_acerbi_population( poplist::PopList, g::Int64, new_id::Vector{Int64},  N::Int64, mu::Float64, 
      conformist_prob::Float64, anti_conformist_prob::Float64, toplist_size::Int64, bottomlist_size::Int64, 
      acerbi_bottomlist::Bool )
Generates a new population using the Acerbi method of conformity and anti-conformity.
"""
function new_acerbi_population( poplist::PopList, g::Int64, new_id::Vector{Int64},  N::Int64, mu::Float64, 
      conformist_prob::Float64, anti_conformist_prob::Float64, toplist_size::Int64, bottomlist_size::Int64, 
      acerbi_bottomlist::Bool )
  toplist = topKlist( poplist[g-1], toplist_size )
  if acerbi_bottomlist
    bottomlist = bottomKlist( poplist[g-1], bottomlist_size )
  end
  result = zeros(Int64,N)
  for i = 1:N
    rnd = rand()
    if rnd < mu  # mutate
      result[i] = new_id[1]
      new_id[1] += 1
    elseif rnd < mu+conformist_prob  # conformist copy
      if poplist[g-1][i] in toplist
        result[i] = poplist[g-1][i]  # copy if in toplist
      else
        result[i] = poplist[g-1][rand(1:N)]
      end
    elseif rnd < mu+conformist_prob+anti_conformist_prob
      if acerbi_bottomlist
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
  return result
end

# TODO:  Check if Kandler's model of conformity agrees with Acerbi's.
@doc """ function acerbi_mixed_conformist_poplist( N::Int64, N_mu::Float64, ngens::Int64, 
    conformist_prob::Float64, anti_conformist_prob::Float64, acerbi_flag::Bool=true;
    toplist_size::Int64=1, bottomlist_size::Int64=1, acerbi_bottomlist::Bool=true,
    burn_in::Float64=2.0, uniform_start::Bool=false, combine::Bool=true )
Return either a list of ngens populations, or a combined population depending on whether "combine" is true.
If acerbi_flag==true, the population will be generated by a mix of random copying, Acerbi conformist 
    copying, and bottomlist or Acerbi anti-conformist copying.
If acerbi_flag==false, the population will be generated by a mix of random, conformist, and
    anti-conformist copying where conformist copying chooses from the toplist, anti-conformist copying
    copies from the bottomlist.  This is the Mesoudi-Lycett method.
bottom==true means use a bottomlist for anti-conformism (only relevant if acerbi_flag == true)
bottom==false means use Acerbi's method of anti-conformism (only relevant if acerbi_flag == true)
n = sample size
N = popsize
N_mu = mutation rate per population, i. e.,  mu*N.  Thus, mu = N_mu/N.
acerbi_flag = true means use the Acerbi/Bentley definition of conformity/anti-conformity
acerbi_flag = false means use the Mesoudi/Lycett definition of conformity/anti-conformity
"""
function acerbi_mixed_conformist_poplist( N::Int64, N_mu::Float64, ngens::Int64, 
    conformist_prob::Float64, anti_conformist_prob::Float64; acerbi_flag::Bool=true,
    toplist_size::Int64=1, bottomlist_size::Int64=1, acerbi_bottomlist::Bool=true,
    burn_in::Float64=2.0, uniform_start::Bool=false, combine::Bool=true )
  if conformist_prob + anti_conformist_prob > 1.0
    error("conformist_prob + anti_conformist_prob must be less than 1.0 in power_mixed_conformist.")
  end
  int_burn_in = Int(round(N*burn_in))
  new_id = zeros(Int64,1)   # new_id[1] is used to assign a new inteter to each innovation
  mu = N_mu/N
  if uniform_start  # All allele values start with the same value.  Start with a selective sweep.
    poplist= Population[ Int64[1 for i = 1:N] ]
    new_id[1] = 2
  else
    poplist= Population[ collect(1:N) ]
    new_id[1] = N+1
  end
  if combine
    pop_result = Population()
  end
  for g = 2:(ngens+int_burn_in)
    result = new_acerbi_population( poplist, g, new_id, N, mu, conformist_prob, anti_conformist_prob, 
        toplist_size, bottomlist_size, acerbi_bottomlist )
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
    #w_conf = fill(1.0,N)
    #w_anti_conf = fill(1.0,N)
    fitness_table = Dict{Int64,Float64}()
    dfe = dfe_neutral
    rp = randperm(N)
    if conformist_prob > 0.0
      conf_fitness_table = freq_scaled_fitness( poplist[g-1], conformist_power, fitness_table )
      new_conf_pop = propsel( poplist[g-1], dfe, conf_fitness_table )[rp]
    end
    if anti_conformist_prob > 0.0
      anti_conf_fitness_table = freq_scaled_fitness( poplist[g-1], anti_conformist_power, fitness_table )
      new_anti_conf_pop = propsel( poplist[g-1], dfe, anti_conf_fitness_table )[rp]
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

function nearly_neutral_power_mixed_conformist_poplist( N::Int64, N_mu::Float64, ngens::Int64,
    conformist_prob::Float64, anti_conformist_prob::Float64; dfe::Function=dfe_neutral,
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
  fitness_table = Dict{Int64,Float64}()  # Fitness table contains the dfe fitness of each id.
  if uniform_start  # All allele values start with the same value.  Start with a selective sweep.
    poplist= Population[ Int64[1 for i = 1:N] ]
    dfe_fitness(1, dfe, fitness_table )  # set fitness of i to be dfe(i).
    new_id = 2
  else
    poplist= Population[ collect(1:N) ]
    for i = 1:N
      # Note that dfe_advantageous(i), dfe_deleterious(i), dfe_mixed(i) do not depend on i, but dfe_mod(i) does depend on i
      dfefit=dfe_fitness(i, dfe, fitness_table )  # set fitness of i to be dfe(i).
    end
    new_id = N+1
  end
  if combine
    pop_result = Population()
  end
  for g = 2:(ngens+int_burn_in)
    rp = randperm(N)
    if conformist_prob > 0.0
      conf_fitness_table = freq_scaled_fitness( poplist[g-1], conformist_power, fitness_table )
      new_conf_pop = propsel( poplist[g-1], dfe, conf_fitness_table )[rp]
    end
    if anti_conformist_prob > 0.0
      anti_conf_fitness_table = freq_scaled_fitness( poplist[g-1], anti_conformist_power, fitness_table )
      new_anti_conf_pop = propsel( poplist[g-1], dfe, anti_conf_fitness_table )[rp]
    end
    result = zeros(Int64,N)
    mutate_count = 0
    conf_count = 0
    anti_conf_count = 0
    random_count = 0
    for i = 1:N
      rnd = rand()
      if rnd < mu  # mutate
        mutate_count += 1
      elseif rnd < mu+conformist_prob  # conformist copy (without replacement)
        conf_count += 1
      elseif rnd < mu+conformist_prob+anti_conformist_prob
        anti_conf_count += 1
      else     # random copy
        random_count += 1
        result[i] = poplist[g-1][rand(1:N)]
      end
    end
    for i = 1:mutate_count
      result[i] = new_id
      dfe_fitness( new_id, dfe, fitness_table )  # Set fitness of new_id
      #println("new_id: ",new_id," fit: ",fitness_table[new_id] )
      new_id += 1
    end
    if conf_count > 0
      conf_fitness_table = freq_scaled_fitness( poplist[g-1], conformist_power, fitness_table )
      result[mutate_count+1:mutate_count+conf_count] = propsel( poplist[g-1], conf_count, dfe, conf_fitness_table )
    end
    if anti_conf_count > 0
      anti_conf_fitness_table = freq_scaled_fitness( poplist[g-1], anti_conformist_power, fitness_table )
      result[mutate_count+conf_count+1:mutate_count+conf_count+anti_conf_count] = propsel( poplist[g-1], anti_conf_count, dfe, anti_conf_fitness_table )
    end
    if random_count > 0
      result[mutate_count+conf_count+anti_conf_count+1:end] = propsel( poplist[g-1], random_count, dfe, fitness_table )
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

