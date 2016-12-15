#using DataStructures
export freq_scaled_fitness, propsel, propsel!

@doc """ function freq_scaled_fitness( pop::Population, c::Float64 )

Return a fitness vector that reflects power conformist or anti-conformist selection on p.
Fitnesses are rescaled by multiplying by the frequency of the corresponding allele to the power c.
The parameter c controls the degree of conformity.
c = 0.0    gives a neutral fitness (the fitness vector w is returned)
c > 0.0    gives a conformist fitness where common elements are favored
c < 0.0    gives an anti-conformist fitness where less elements are favored
"""
#=
function freq_scaled_fitness( pop::Population, w::Vector{Float64}, c::Float64 )
  frequency = DataStructures.counter(Int64)
  for p in pop 
    push!( frequency, p )
  end
  n = length(pop)
  [ w[pop[i]]*frequency[pop[i]]^c for i in 1:n ]
end
=#
@doc """
Same as previous version except that fitness is stored in a dictionary
"""
function freq_scaled_fitness( pop::Population, c::Float64, fitness_table::Dict{Int64,Float64} )
  frequency = DataStructures.counter(Int64)
  for p in pop 
    push!( frequency, p )
  end
  n = length(pop)
  new_w = Dict{Int64,Float64}()
  for i in 1:n
    if get(fitness_table,pop[i],-1.0) == -1.0
      new_w[pop[i]] = frequency[pop[i]]^c
    else 
      new_w[pop[i]] = fitness_table[pop[i]]*frequency[pop[i]]^c
    end
  end
  new_w
end

#=
@doc """ function propsel( pop::Population, fitness::Vector{Float64}  )
Apply proportional selection to Population pop using fitness, 
and return the result.  
"""
function propsel( pop::Population, fitness::Vector{Float64}  )
  new_pop = deepcopy(pop)
  propsel!( new_pop, fitness )
  new_pop
end
=#

@doc """ function propsel( pop::Population, fitness::Function )
Apply proportional selection to Population pop using fitness, 
and return the result.  
"""
function propsel( pop::Population, dfe::Function, fitness_table::Dict{Int64,Float64} )
  new_pop = deepcopy(pop)
  propsel!( new_pop, dfe, fitness_table )
  new_pop
end

@doc """ function propsel( pop::Population, fitness::Function )
Apply proportional selection to Population pop using fitness to generate
and return a new selected population of size n.
"""
function propsel( pop::Population, n::Int64, dfe::Function, fitness_table::Dict{Int64,Float64} )
  new_pop = zeros(Int64,n)
  fmax = 0.0
  for p in pop
    fitp = dfe_fitness( p, dfe, fitness_table )
    if fitp > fmax
      fmax = fitp
    end
  end
  if fmax == 0.0
    # all elements have fitness zero
    return
  end

  N = length(pop)
  k = 0
  while k < n
    i = rand(1:N)
    w = dfe_fitness(pop[i],dfe, fitness_table ) / fmax
    if rand() < w
      new_pop[k + 1] = pop[i]
      k += 1
    end
  end
  new_pop
end

#=
@doc """function propsel!(p::Population, fitness::Vector{Float64} )
Conduct proportional selection in-place.
"""
function propsel!( pop::Population, fitness::Vector{Float64}  )
  fmax = maximum(fitness)
  if fmax == 0
    # all elements have fitness zero
    return
  end

  n = length(pop)
  selected = zeros(Int64, n)
  k = 0
  while k < n
    i = rand(1:n)
    w = fitness[pop[i]] / fmax
    if rand() < w
      selected[k + 1] = i
      k += 1
    end
  end
  pop[:] = [ pop[selected[i]] for i = 1:n ]
end
=#

@doc """function propsel!(p::Population, fitness::Function, dfe::Function )
Conduct proportional selection in-place.
"""
function propsel!( pop::Population, dfe::Function, fitness_table::Dict{Int64,Float64} )
  fmax = 0.0
  for p in pop
    fitp = dfe_fitness( p, dfe, fitness_table )
    if fitp > fmax
      fmax = fitp
    end
  end
  if fmax == 0.0
    # all elements have fitness zero
    return
  end

  n = length(pop)
  selected = zeros(Int64, n)
  k = 0
  while k < n
    i = rand(1:n)
    w = dfe_fitness(pop[i], dfe, fitness_table ) / fmax
    if rand() < w
      selected[k + 1] = i
      k += 1
    end
  end
  pop[:] = [ pop[selected[i]] for i = 1:n ]
end
