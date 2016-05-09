using DataStructures
typealias Population Array{Int64,1}
typealias PopList Array{Population,1}

function power_confirmist_poplist( N::Int64, mu::Float64, ngens::Int64; burn_in::Int64=0, conformist_power::Float64=0.0 )
  poplist= Population[ collect(1:N) ]
  new_id = N+1
  for g = 2:(ngens+burn_in)
    fitness = freq_scaled_fitness( poplist[g-1], conformist_power )
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

function acerbi_conformist_poplist( N::Int64, mu::Float64, ngens::Int64, K::Int64, C::Float64; burn_in::Int64=0 )
  poplist= Population[ collect(1:N) ]
  new_id = N+1
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

function acerbi_anti_conformist_poplist( N::Int64, mu::Float64, ngens::Int64, K::Int64, C::Float64; burn_in::Int64=0 )
  poplist= Population[ collect(1:N) ]
  new_id = N+1
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


