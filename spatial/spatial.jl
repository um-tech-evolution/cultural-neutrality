# Spatial structure simulation with horizontal transfer

using Distributions
typealias Population Array{Int64,1}
typealias PopList Array{Population,1}

type variant_type
  parent::Int64   # The variant that gave rise to this variant
  innovation::Int64   # The innovation that gave rise to this variant
  fitness::Float64    # The fitness of this variant
end

type subpop_type
  fitness_inc::Float64
end

empty_variant = variant_type(-1,-1,0.0)
#vtbl = Dict{Int64,variant_type}()

@doc """ function wf_simmulation()
  Wright-Fisher model simulation (as opposed to Moran model)
  Parameters:
    N     MetaPopulation size
    m     number of subpopulations   # for now, subpopulation size = N/m
    mu    innovation probability
    copy_err_prob    copy error probability
    ngens number of generations
    copy_dfe   Distribution of selection coefficients for the reduction in fitness during copy
    innov_dfe  Distribution of selection coefficients for innovations
    horiz_dfe  Distribution of selection coefficients for the reduction in fitness during horiz transfer
    variant_table Keeps track fitnesses and variant parent and innovation ancestor
"""
function spatial_simulation( N::Int64, m::Int64, mu::Float64, copy_err_prob::Float64, ngens::Int64, 
    num_emmigrants::Int64, copy_dfe::Function, innov_dfe::Function, variant_table::Dict{Int64,variant_type};
    forward::Bool=true, neg_select::Bool=true )
  id = Int[0]
  n = Int(floor(N/m))    # size of subpopulations
  println("n: ",n)
  # Initialize
  pop_list = Vector{PopList}()
  subpops = PopList()
  for j = 1:m
    Base.push!( subpops, Population() )
    for i = 1:n
      Base.push!( subpops[j], new_innovation( id, innov_dfe, variant_table ) )
    end
    println("subpops[",j,"]: ",subpops[j] )
  end
  Base.push!(pop_list,deepcopy(subpops))
  for g = 2:ngens
    #println("g: ",g)
    for j = 1:m
      for i = 1:n
        cp = copy_parent( pop_list[g-1][j][i], copy_err_prob, id, copy_dfe, variant_table )
        #println("j: ",j,"  i: ",i,"  pl: ",pop_list[g-1][j][i],"  cp: ",cp)
        subpops[j][i] = cp
      end
      subpops[j] = propsel( subpops[j], n, variant_table )
      println("subpops[",j,"]: ",map(x->variant_table[x].fitness,subpops[j] ))
      horiz_transfer_circular!( N, m, num_emmigrants, subpops, variant_table )
    end
    Base.push!(pop_list,deepcopy(subpops))
  end
  pop_list
end

function fitness( v::Int64, innov::Int64, dfe::Function, variant_table::Dict{Int64,variant_type} )
  vt = get( variant_table, v, -1 )
  if vt == -1
    fit = dfe() 
    vt = variant_type( v, innov, fit )
    variant_table[v] = vt
    return fit
  else
    return vt.fitness
  end
end

function new_innovation( id::Vector{Int64}, dfe::Function, variant_table::Dict{Int64,variant_type} )
  #println("new innovation  dfe: ",dfe)
  i = id[1]
  fit = fitness( i, i, dfe, variant_table )
  #println("i: ",i,"  fit: ",fit)
  id[1] += 1
  i
end

@doc """  copy_parent()
  Note that the function cdfe produces an incremental change in fitness rather than a new fitness.
"""
function copy_parent( v::Int64, copy_err_prob::Float64, id::Vector{Int64}, cdfe::Function, variant_table::Dict{Int64,variant_type} )
  if rand() < copy_err_prob
    i = id[1]
    vt = variant_table[v]
    ffit = cdfe()
    new_fit = vt.fitness + ffit
    #println("copy_parent: ffit: ",ffit,"  new_fit: ",new_fit)
    variant_table[i] = variant_type(v,vt.innovation,new_fit)
    id[1] += 1
    return i
  else
    return v
  end
end  

@doc """ horiz_transfer_circular!()
  Transfers variants between subpopulations in a circular fashion (either forward or backward).
  Elements to be transfered are selected by proportional selection.
  Elements to be replaced can be random or selected by reverse proportional selection depending on the flag neg_select.
  subpops is modified by this function (as a side effect)
"""
function horiz_transfer_circular!( N::Int64, m::Int64, num_emmigrants::Int64, subpops::PopList, variant_table::Dict{Int64,variant_type};
     forward::Bool=true, neg_select::Bool=true )
  n = Int(floor(N/m))    # size of subpopulations
  emmigrants = PopList()
  for j = 1:m
    Base.push!( emmigrants, propsel( subpops[j], num_emmigrants, variant_table ) )
  end
  for j = 1:m
    #println("j: ",j,"  j%m+1: ",j%m+1,"  (j+m-2)%m+1: ",(j+m-2)%m+1)
    if neg_select  # use reverse proportional selection to delete elements by negative fitness
      pop_after_deletion = reverse_propsel(subpops[j],num_emmigrants,variant_table)
    else  # delete random elements to delete
      pop_after_deletion = subpops[j][1:(n-num_emmigrants)]
    end
    if forward
      subpops[j] = append!( pop_after_deletion, emmigrants[(j+m-2)%m+1] )
    else
      subpops[j] = append!( pop_after_deletion, emmigrants[j%m+1] )
    end
  end
  emmigrants  # perhaps should be the modified subpops
end

include("../src/NeutralCulturalEvolution.jl")
include("../spatial/propsel.jl")
if length(ARGS) == 0
  simname = "configs/example1"
else
  simname = ARGS[1]
end
include("$(simname).jl")
println("simname: ",simname)
#stream = open("$(simname).csv","w")
#println("stream: ",stream)


if isdefined(:simtype)
  pop_list = spatial_simulation( N, m, mu, copy_err_prob, ngens, ne, cdfe, idfe, vtbl )
  for i = 1:ngens
    println(pop_list[i])
  end
end

  
function init()
  global vtbl = Dict{Int64,variant_type}()
  global idfe() = 1.0+rand(Distributions.Gamma(1.0,0.1))
  global cdfe() = -rand(Distributions.Gamma(0.2,0.001))
  global N = 20
  global m = 4
  global mu = 0.2
  global ne = 3  # num_emmigrants 
  global ngens = 4
end  
# Sample call to main function
# spatial_simulation( N, m, mu, ngens, ne, cdfe, idfe, vtbl )
