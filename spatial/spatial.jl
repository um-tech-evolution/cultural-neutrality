# Spatial structure simulation with horizontal transfer
export spatial_simulation, fitness

#=  has been moved to types.jl
using Distributions
typealias Population Array{Int64,1}
typealias PopList Array{Population,1}

type variant_type
  parent::Int64   # The variant that gave rise to this variant
  innovation::Int64   # The innovation that gave rise to this variant
  fitness::Float64    # The fitness of this variant
  subpop_index::Int64  # index of the containing subpopulation
  attributes::Vector{Float64}   # attributes of the variant
end

type subpop_type
  simtype::Int64
  fitness_inc::Float64
  ideal::Vector{Float64}   # Ideal values for attributes in the environment of this subpop
end
=#

empty_variant = variant_type(-1,-1,0.0,0,Vector{Float64}())
#vtbl = Dict{Int64,variant_type}()

@doc """ function spatial_simulation()
  Wright-Fisher model simulation (as opposed to Moran model)
  Parameters:
    N     MetaPopulation size
    m     number of subpopulations   # for now, subpopulation size = N/m
    mu    innovation probability
    copy_err_prob    copy error probability
    ngens number of generations
    num_emmigrants   number of emmigrants in horizontal transfer
    num_attributes   number of quantitative attributes of a variant
    copy_dfe   Distribution of selection coefficients for the reduction in fitness during copy
    innov_dfe  Distribution of selection coefficients for innovations
    horiz_dfe  Distribution of selection coefficients for the reduction in fitness during horiz transfer
    variant_table Keeps track fitnesses and variant parent and innovation ancestor
    quantitative==true means individuals have quantitative attributes, fitness computed by distance from ideal
    forward==true  means that horizontal transfer is done in a forward circular fashion
    neg_select==true  means that reverse proportional selection is used to select individuals to delete in horiz trans
"""
function spatial_simulation( N::Int64, num_subpops::Int64, mu::Float64, copy_err_prob::Float64, ngens::Int64, burn_in::Float64,
    num_emmigrants::Int64, num_attributes::Int64, normal_stddev::Float64, copy_dfe::Function, innov_dfe::Function; 
    #variant_table::Dict{Int64,variant_type};
    quantitative::Bool=true, forward::Bool=true, neg_select::Bool=true )
  variant_table = Dict{Int64,variant_type}()
  println("spatial simulation: quantitative: ",quantitative)
  subpop_properties = initialize_subpop_properties(num_subpops,num_attributes)
  #println("subpop_properties: ",subpop_properties)
  #int_burn_in = Int(round(burn_in*N+50.0))
  int_burn_in = Int(round(burn_in*N))   # reduce for testing
  id = Int[1]
  n = Int(floor(N/num_subpops))    # size of subpopulations
  println("N: ",N,"  num_subpops: ",num_subpops,"  n: ",n,"  num_attributes: ",num_attributes,"  ngens: ",ngens)
  # Initialize
  cumm_means = zeros(Float64,num_subpops)
  cumm_variances = zeros(Float64,num_subpops)
  cumm_attr_vars = zeros(Float64,num_subpops)
  pop_list = Vector{PopList}()
  subpops = PopList()
  for j = 1:num_subpops
    Base.push!( subpops, Population() )
    for i = 1:n
      Base.push!( subpops[j], new_innovation( id, innov_dfe, j, num_attributes, variant_table, subpop_properties ) )
    end
    #println("subpops[",j,"]: ",subpops[j] )
  end
  Base.push!(pop_list,deepcopy(subpops))
  for g = 2:ngens+int_burn_in
    #println("g: ",g)
    for j = 1:num_subpops
      for i = 1:n
        cp = copy_parent( pop_list[g-1][j][i], copy_err_prob, id, j, copy_dfe, normal_stddev, variant_table, subpop_properties )
        #println("j: ",j,"  i: ",i,"  pl: ",pop_list[g-1][j][i],"  cp: ",cp)
        subpops[j][i] = cp
      end
      subpops[j] = propsel( subpops[j], n, variant_table )
      #println("subpops[",j,"]: ",map(x->variant_table[x].fitness,subpops[j] ))
      horiz_transfer_circular!( N, num_subpops, num_emmigrants, subpops, variant_table )
    end
    Base.push!(pop_list,deepcopy(subpops))
    #print_pop(STDOUT,subpops,variant_table)
    if g > int_burn_in
      (mmeans,vvars) = means(subpops,variant_table)
      cumm_means += mmeans
      cumm_variances += vvars
      avars = attr_vars(subpops,variant_table)
      cumm_attr_vars += avars
      #println("cumm_means: ",cumm_means)
      #println("cumm_variances: ",cumm_variances)
    end
  end
  cumm_means /= ngens
  cumm_variances /= ngens
  cumm_attr_vars /= ngens
  println("cumm_means: ",cumm_means)
  println("cumm_variances: ",cumm_variances)
  println("cumm_attr_vars: ",cumm_attr_vars)
  #=
  println("pop_list: ")
  for i = int_burn_in+1:ngens+int_burn_in
    println(pop_list[i])
  end
  =#
  return mean(cumm_means), mean(cumm_variances), mean(cumm_attr_vars)
end

function fitness( v::Int64, innov::Int64, dfe::Function, subpop_index::Int64, 
    variant_table::Dict{Int64,variant_type} )
  #println("fitness: v: ",v)
  vt = get( variant_table, v, -1 )
  if vt == -1
    fit = dfe() 
    vt = empty_variant
    vt.parent = v
    vt.innovation = innov
    vt.fitness = fit
    vt.subpop_index = subpop_index
    #vt = variant_type( v, innov, fit, variant_table[v].subpop_index, variant_table[v].attributes )
    variant_table[v] = vt
    return fit
  else
    return vt.fitness
  end
end

function fitness( attributes::Vector{Float64}, ideal::Vector{Float64} )
  if length(attributes) != length(ideal)
    error("length(attributes) must equal length(ideal) in fitness")
  end
  sum = 0.0
  for k = 1:length(attributes)
    sum += abs( attributes[k] - ideal[k] )
  end
  return 1.0-sum/length(attributes)
end

function new_innovation( id::Vector{Int64}, dfe::Function, subpop_index::Int64, num_attributes::Int64,
    variant_table::Dict{Int64,variant_type}, subpop_properties; quantitative::Bool=true )
  i = id[1]
  if quantitative
    #println("new innovation i: ",i,"  subpop_index:",subpop_index,"  num_attributes: ",num_attributes )
    variant_table[i] = variant_type( i, i, 0.0, subpop_index, rand(num_attributes) )
    #println("variant_table[i]: ",variant_table[i])
    variant_table[i].fitness = fitness( variant_table[i].attributes, subpop_properties[subpop_index].ideal )  
  else
    #println("new innovation i: ",i,"  dfe: ",dfe)
    fit = fitness( i, i, dfe, subpop_index, variant_table )
    #println("i: ",i,"  fit: ",fit)
  end
  id[1] += 1
  i
end

@doc """  copy_parent()
  Note that the function cdfe produces an incremental change in fitness rather than a new fitness.
"""
function copy_parent( v::Int64, copy_err_prob::Float64, id::Vector{Int64}, subpop_index::Int64, cdfe::Function, 
    normal_stddev::Float64, variant_table::Dict{Int64,SpatialEvolution.variant_type},
    subpop_properties::Vector{SpatialEvolution.subpop_type}; quantitative::Bool=true )
  if quantitative   # temporary fix
    copy_err_prob = 1.0
  end
  if rand() < copy_err_prob
    i = id[1]
    vt = variant_table[v]
    if quantitative  # spatial structure by deviation of attributes from ideal
      vt.attributes = mutate_attributes( vt.attributes, normal_stddev )
      new_fit = fitness( vt.attributes, subpop_properties[subpop_index].ideal )
      #=
      distance = 0.0
      println("length(vt.attributes): ",length(vt.attributes),"  length(ideal): ",length(subpop_properties[subpop_index].ideal))
      for i = 1:length(vt.attributes)
        distance += abs(vt.attributes[i] - subpop_properties[subpop_index].ideal[i])
      end
      new_fit = 1.0-distance/length(vt.attributes)
      =#
    else   # spatial structure by fitness increment
      ffit = cdfe()
      new_fit = vt.fitness + ffit
    end
    #println("copy_parent i: ",i,"  quantitative: ",quantitative,"  new_fit: ",new_fit)
    variant_table[i] = variant_type(v,vt.innovation,new_fit,vt.subpop_index,vt.attributes)
    #println("v: ",v,"  i: ",i,"  new_fit: ",new_fit,"  vtbl[i]: ",variant_table[i].fitness)
    id[1] += 1
    return i
  else
    return v
  end
end  

function mutate_attributes( attributes::Vector{Float64}, normal_stddev::Float64 )
  #stddev = normal_stddev()   # Standard deviation of mutational perturbations
  #println("mutate attributes  normal_stddev: ",normal_stddev)
  attributes = min(1.0,max(0.0,attributes+normal_stddev*randn(length(attributes))))
  #println("attributes: ",attributes)
  return attributes
end

function initialize_subpop_properties( num_subpops::Int64, num_attributes::Int64) 
    #subpop_properties::Vector{subpop_type} )
  fit_inc = 0.1
  subpop_properties = [ subpop_type( fit_inc, zeros(Float64,num_attributes) ) for j = 1:num_subpops ]
  min_value = 0.45
  max_value = 0.55
  for j = 1:num_subpops
    for k = 1:num_attributes
      if min_value != max_value
        subpop_properties[j].ideal[k] = min_value+rand()*(max_value-min_value)
      else
        subpop_properties[j].ideal[k] = min_value
      end
    end
  end
  return subpop_properties
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

function print_subpop( subpop::Vector{Int64}, variant_table::Dict{Int64,variant_type} )
  "[ $(join([ @sprintf(" %5d:%5.4f",vt,variant_table[vt].fitness) for vt in subpop ]))]"
end

function print_pop( stream::IO, subpops::PopList, variant_table::Dict{Int64,variant_type} )
  for sp in subpops
    print(stream,print_subpop(sp,variant_table))
  end
  println(stream)
end

using DataFrames
# compute and save statistics about subpopulations and populations

function means( subpops::PopList, variant_table::Dict{Int64,variant_type} )
  fit(v) = variant_table[v].fitness
  means = [ mean(map(fit,s)) for s in subpops ]
  vars = [ var(map(fit,s)) for s in subpops ]
  return means, vars
end

function attr_vars( subpops::PopList, variant_table::Dict{Int64,variant_type} )
  num_attributes = length(variant_table[1].attributes)
  #=
  println("attr_vars: num_attributes: ",num_attributes)
  for s in subpops
    println(s," fitness: ",[variant_table[v].fitness for v in s ],"  variance: ",var([variant_table[v].fitness for v in s ]))
    for j = 1:num_attributes
      println(j,"  attribute: ",[variant_table[v].attributes[j] for v in s],"  variance: ",var([variant_table[v].attributes[j] for v in s]) )
    end
  end
  =#
  ave_vars = zeros(Float64,length(subpops))
  i = 1
  for s in subpops
    att_vars = [ var([variant_table[v].attributes[j] for v in s]) for j =1:num_attributes]
    #println(s," att_vars: ",att_vars)
    ave_vars[i] = mean(att_vars)
    i += 1
  end
  #println("ave_vars: ",ave_vars)
  return ave_vars
end
 
function init()
  global vtbl = Dict{Int64,variant_type}()
  global subpop_properties = subpop_type[ subpop_type(0.1,zeros(Float64,m)) ]
  global idfe() = 1.0+rand(Distributions.Gamma(1.0,0.1))
  global cdfe() = -rand(Distributions.Gamma(0.2,0.001))
  global N = 20
  global m = 4
  global mu = 0.2
  global cperr = 0.2
  global ne = 3  # num_emmigrants 
  global ngens = 4
end  
# Sample call to main function
# spatial_simulation( N, m, mu, cperr, ngens, ne, cdfe, idfe, vtbl, subpop_properties )
