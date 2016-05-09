using DataStructures
typealias Population Array{Int64,1}
typealias PopList Array{Array{Int64,1},1}

@doc """ function freq_scaled_fitness( pop::Population, pwr::Float64 )

Return a fitness vector that reflects conformist or anti-conformist selection on pop.
The paramter pwr controls the degree of conformity.
pwr = 0.0    gives a neutral fitness
pwr > 0.0    gives a conformist fitness where common elements are favored
pwr < 0.0    gives an anti-conformist fitness where less elements are favored
"""

function freq_scaled_fitness( pop::Population, pwr::Float64 )
  cntr = counter(Int64)
  for i in pop 
    push!( cntr, i )
  end
  n = length(pop)
  [ cntr[i]^pwr/n for i in pop ]
end


#=
@doc """propsel(p::Population, fitness::Vector{Float64} )

Returns proportional selection applied to Population  pop
"""
function propsel( pop::Population, fitness::Vector{Float64}  )
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
    w = fitness[i] / fmax
    if rand() < w
      selected[k + 1] = i
      k += 1
    end
  end
  [ pop[selected[i]] for i = 1:n ]
end
=#

@doc """ function propsel( pop::Population, fitness::Vector{Float64}  )

Apply proportional selection to Population pop using fitness, 
and return the result.  
"""

function propsel( pop::Population, fitness::Vector{Float64}  )
  new_pop = deepcopy(pop)
  propsel!( new_pop, fitness )
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
    w = fitness[i] / fmax
    if rand() < w
      selected[k + 1] = i
      k += 1
    end
  end
  pop[:] = [ pop[selected[i]] for i = 1:n ]
end
