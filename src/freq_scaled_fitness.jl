using DataStructures
export freq_scaled_fitness, propsel, propsel!

@doc """ function freq_scaled_fitness( pop::Population, c::Float64 )

Return a fitness vector that reflects conformist or anti-conformist selection on p.
Fitnesses are rescaled by multiplying by the frequency of the corresponding allele to the power c.
The w parameter gives the fitnesses of the population.  
The parameter c controls the degree of conformity.
c = 0.0    gives a neutral fitness (the fitness vector w is returned)
c > 0.0    gives a conformist fitness where common elements are favored
c < 0.0    gives an anti-conformist fitness where less elements are favored
"""

function freq_scaled_fitness( pop::Population, w::Vector{Float64}, c::Float64 )
  frequency = counter(Int64)
  for p in pop 
    push!( frequency, p )
  end
  n = length(pop)
  [ w[i]*frequency[pop[i]]^c for i in 1:n ]
end

@doc """ function propsel( pop::Population, fitness::Vector{Float64}  )
Apply proportional selection to Population pop using fitness, 
and return the result.  
"""
function propsel( pop::Population, fitness::Vector{Float64}  )
  new_pop = deepcopy(pop)
  propsel!( new_pop, fitness )
  new_pop
end

@doc """ function propsel( pop::Population, fitness::Function )
Apply proportional selection to Population pop using fitness, 
and return the result.  
"""
function propsel( pop::Population, fitness::Function, dfe::Function )
  new_pop = deepcopy(pop)
  propsel!( new_pop, fitness, dfe )
  new_pop
end

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

@doc """function propsel!(p::Population, fitness::Function, dfe::Function )
Conduct proportional selection in-place.
"""
function propsel!( pop::Population, fitness::Function, dfe::Function )
  fmax = 0.0
  for p in pop
    fitp = fitness( p, dfe )
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
    w = fitness(pop[i],dfe) / fmax
    if rand() < w
      selected[k + 1] = i
      k += 1
    end
  end
  pop[:] = [ pop[selected[i]] for i = 1:n ]
end
