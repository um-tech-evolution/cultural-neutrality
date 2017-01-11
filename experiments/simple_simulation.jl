# simple simulation to illustrate selection intensity

include("../src/freq_scaled_fitness.jl")

function dfe_simple( x:: Int64, N::Int64; fit_inc::Float64=1.0/N )
  if x == 1 
    return 1.0 + fit_inc 
  else
    return 1.0
  end
end


function s_poplist( N::Int64, ngens::Int64, k::Int64, dfe::Function )
  global fitness_table = Dict{Int64,Float64}()
  new_pop = vcat([1 for i = 1:k],[ 2 for i = 1:(N-k)])
  poplist= Population[ new_pop ]
  fit = dfe_fitness(1, dfe, fitness_table )  # set fitness of 1 to be dfe(1).
  fit = dfe_fitness(2, dfe, fitness_table )  # set fitness of 1 to be dfe(1).
  #println("new pop: ",new_pop)
  for g = 2:ngens
    new_pop = propsel( poplist[g-1], dfe, fitness_table )
    #println("new pop: ",new_pop)
    Base.push!(poplist,new_pop)
    u = unique(new_pop)
    #println("unique: ",u)
    if length(u) == 1
      return (u[1],g)
    end
  end
  poplist
end

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

function run_trials_rel( N::Int64, k::Int64, T::Int64, fit_inc::Float64=0.0 )
  ngens = 100*N
  function dfe_s(x)
    if x == 1
      return 1.0 + fit_inc/N 
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
    pl = s_poplist(N,ngens,k,dfe_s)
    #println("pl: ",pl)
    count_1 += pl[1]==1?1:0
    count_2 += pl[1]==2?1:0
    gens_1 += pl[1]==1?pl[2]:0
    gens_2 += pl[1]==2?pl[2]:0
  end
  FT = Float64(T)
  println("counts: 1: ",count_1/FT,"  2: ",count_2/FT)
  println("gens:   1: ", gens_1/FT,"  2: ",gens_2/FT)
end

function run_trials_fixed( N::Int64, k::Int64, T::Int64, s::Float64=0.0 )
  ngens = 100*N
  function dfe_s(x)
    if x == 1
      return  s
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
    pl = s_poplist(N,ngens,k,dfe_s)
    #println("pl: ",pl)
    count_1 += pl[1]==1?1:0
    count_2 += pl[1]==2?1:0
    gens_1 += pl[1]==1?pl[2]:0
    gens_2 += pl[1]==2?pl[2]:0
  end
  FT = Float64(T)
  println("counts: 1: ",count_1/FT,"  2: ",count_2/FT)
  println("gens:   1: ", gens_1/FT,"  2: ",gens_2/FT)
end
