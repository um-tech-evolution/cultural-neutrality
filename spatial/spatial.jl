# Spatial structure simulation with horizontal transfer
export spatial_simulation, fitness

empty_variant = variant_type(-1,-1,0.0,0,Vector{Float64}())
#vtbl = Dict{Int64,variant_type}()

@doc """ function spatial_simulation()
  Wright-Fisher model simulation (as opposed to Moran model)
  Parameters:
    N     MetaPopulation size
    m     number of subpopulations   # for now, subpopulation size = N/m
    mu    innovation probability
    ngens number of generations
    num_emmigrants   number of emmigrants in horizontal transfer
    num_attributes   number of quantitative attributes of a variant
    copy_dfe   Distribution of selection coefficients for the reduction in fitness during copy
    innov_dfe  Distribution of selection coefficients for innovations
    horiz_dfe  Distribution of selection coefficients for the reduction in fitness during horiz transfer
    variant_table Keeps track fitnesses and variant parent and innovation ancestor
    quantitative==true means individuals have quantitative attributes, fitness computed by distance from ideal
    forward==true  means that horizontal transfer is done in a forward circular fashion
    extreme==true  means that horizontal transfer is done in a forward circular fashion
    neg_select==true  means that reverse proportional selection is used to select individuals to delete in horiz trans
"""
#=
function spatial_simulation( N::Int64, num_subpops::Int64, mu::Float64, ngens::Int64, burn_in::Float64,
    num_emmigrants::Int64, num_attributes::Int64, normal_stddev::Float64; 
    quantitative::Bool=true, forward::Bool=true, neg_select::Bool=true,  
    horiz_select::Bool=true,
    extreme_variation::Bool=false,      # Ideal values randomly alternate between a high and a low values
    circular_variation::Bool=false )    # Vary ideal values in a circular fashion
=#
function spatial_simulation( sr::SpatialEvolution.spatial_result_type )
  variant_table = Dict{Int64,variant_type}()
  #println("spatial simulation: quantitative: ",quantitative)
  #println("sim circular_variation: ",sr.circular_variation,"  extreme_variation: ",sr.extreme_variation)
  if sr.num_env_subpops == 0
    ideal_properties = initialize_ideal_properties(sr.num_subpops,sr.num_attributes,
      extreme_variation=sr.extreme_variation, circular_variation=sr.circular_variation)
  else
    ideal_properties = initialize_ideal_properties(sr.num_env_subpops,sr.num_attributes,
      extreme_variation=sr.extreme_variation, circular_variation=sr.circular_variation)
  end
  #println("ideal_properties: ",ideal_properties)
  #int_burn_in = Int(round(sr.burn_in*N+50.0))
  int_burn_in = Int(round(sr.burn_in*sr.N))   # reduce for testing
  id = Int[1]
  n = Int(floor(sr.N/sr.num_subpops))    # size of subpopulations
  #println("N: ",sr.N,"  num_subpops: ",sr.num_subpops,"  n: ",n,"  num_attributes: ",sr.num_attributes,"  ngens: ",sr.ngens)
  cumm_means = zeros(Float64,sr.num_subpops)
  cumm_variances = zeros(Float64,sr.num_subpops)
  cumm_attr_vars = zeros(Float64,sr.num_subpops)
  pop_list = Vector{PopList}()
  subpops = PopList()
  for j = 1:sr.num_subpops
    Base.push!( subpops, Population() )
    for i = 1:n
      Base.push!( subpops[j], new_innovation( id, 
          j, sr.num_attributes, variant_table, ideal_properties ) )
    end
    #println("subpops[",j,"]: ",subpops[j] )
  end
  Base.push!(pop_list,deepcopy(subpops))
  for g = 2:sr.ngens+int_burn_in
    #println("g: ",g)
    for j = 1:sr.num_subpops
      for i = 1:n
        cp = copy_parent( pop_list[g-1][j][i], id, j, sr.mu, 
          sr.normal_stddev, variant_table, ideal_properties )
        #println("j: ",j,"  i: ",i,"  pl: ",pop_list[g-1][j][i],"  cp: ",cp)
        subpops[j][i] = cp
      end
      subpops[j] = propsel( subpops[j], n, variant_table )
      #println("subpops[",j,"]: ",map(x->variant_table[x].fitness,subpops[j] ))
      if g%2==0
        horiz_transfer_circular!( sr.N, sr.num_subpops, sr.ne, subpops, variant_table, forward=true, neg_select=sr.horiz_select, emmigrant_select=sr.horiz_select )
      else
        horiz_transfer_circular!( sr.N, sr.num_subpops, sr.ne, subpops, variant_table, forward=false, neg_select=sr.horiz_select, emmigrant_select=sr.horiz_select )
      end
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
  cumm_means /= sr.ngens
  cumm_variances /= sr.ngens
  cumm_attr_vars /= sr.ngens
  #println("fitness mean: ",mean(cumm_means),"  variance: ",mean(cumm_variances),"  attribute_variance: ",mean(cumm_attr_vars))
  sr.fitness_mean = mean(cumm_means)
  sr.fitness_variance = mean(cumm_variances)
  sr.attribute_variance = mean(cumm_attr_vars)
  return sr
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

function new_innovation( id::Vector{Int64}, 
    #dfe::Function, 
    subpop_index::Int64, num_attributes::Int64,
    variant_table::Dict{Int64,variant_type}, ideal_properties; quantitative::Bool=true )
  i = id[1]
  if quantitative
    #println("new innovation i: ",i,"  subpop_index:",subpop_index,"  num_attributes: ",num_attributes )
    variant_table[i] = variant_type( i, i, 0.0, subpop_index, rand(num_attributes) )
    #println("variant_table[i]: ",variant_table[i])
    variant_table[i].fitness = fitness( variant_table[i].attributes, ideal_properties[subpop_index].ideal )  
  end
  id[1] += 1
  i
end

@doc """  copy_parent()
  Note that the function cdfe produces an incremental change in fitness rather than a new fitness.
"""
function copy_parent( v::Int64, id::Vector{Int64}, subpop_index::Int64, mu::Float64,
    #cdfe::Function, 
    normal_stddev::Float64, variant_table::Dict{Int64,SpatialEvolution.variant_type},
    ideal_properties::Vector{SpatialEvolution.ideal_type}; quantitative::Bool=true )
  i = id[1]
  vt = variant_table[v]
  if quantitative  # spatial structure by deviation of attributes from ideal
    vt.attributes = mutate_attributes( vt.attributes, normal_stddev )
    if rand() < mu
      innovate_attribute( vt.attributes, subpop_index, ideal_properties )
    end
    new_fit = fitness( vt.attributes, ideal_properties[subpop_index].ideal )
    #=
    distance = 0.0
    println("length(vt.attributes): ",length(vt.attributes),"  length(ideal): ",length(ideal_properties[subpop_index].ideal))
    for i = 1:length(vt.attributes)
      distance += abs(vt.attributes[i] - ideal_properties[subpop_index].ideal[i])
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
end  

function mutate_attributes( attributes::Vector{Float64}, normal_stddev::Float64 )
  #stddev = normal_stddev()   # Standard deviation of mutational perturbations
  #println("mutate attributes  normal_stddev: ",normal_stddev)
  attributes = min(1.0,max(0.0,attributes+normal_stddev*randn(length(attributes))))
  #println("attributes: ",attributes)
  return attributes
end

function innovate_attribute( attributes::Vector{Float64}, subpop_index::Int64, ideal_properties::Vector{SpatialEvolution.ideal_type} )
  j = rand(1:length(attributes))   # Choose a random attribute
  #println("j: ",j,"  attribute: ",attributes[j],"  ideal: ",ideal_properties[subpop_index].ideal[j])
  attributes[j] += rand()*abs(attributes[j] - ideal_properties[subpop_index].ideal[j])*(ideal_properties[subpop_index].ideal[j]-attributes[j])
  #println("j: ",j,"  attribute: ",attributes[j],"  ideal: ",ideal_properties[subpop_index].ideal[j])
end 

function initialize_ideal_properties( num_subpops::Int64, num_attributes::Int64; circular_variation::Bool=false, extreme_variation::Bool=false ) 
  println("init circular_variation: num_subpops: ",num_subpops,"  circ_var: ",circular_variation,"  extreme_var: ",extreme_variation)
  ideal_properties = [ ideal_type( zeros(Float64,num_attributes) ) for j = 1:num_subpops ]
  if !circular_variation && !extreme_variation  # random variation---no relationship to subpop number
    #println("init circular_variation: ",circular_variation,"  extreme_variation: ",extreme_variation)
    min_value = 0.45
    max_value = 0.55
    for j = 1:num_subpops
      for k = 1:num_attributes
        if min_value != max_value
          ideal_properties[j].ideal[k] = min_value+rand()*(max_value-min_value)
        else
          ideal_properties[j].ideal[k] = min_value
        end
      end
    end
  elseif circular_variation && !extreme_variation
    # TODO:  move these parameters to the configuration file
    low_value = 0.2
    high_value = 0.8
    range = 0.2
    increment = 2.0*(high_value-low_value)/num_subpops
    mid = Int(floor(num_subpops/2))
    for j = 1:mid
      for k = 1:num_attributes
        ideal_properties[j].ideal[k] = low_value+increment*(j-1)+(rand()*range-0.5*range)
      end
    end
    for j = (mid+1):num_subpops
      for k = 1:num_attributes
        ideal_properties[j].ideal[k] = high_value-increment*(j-mid-1)+(rand()*range-0.5*range)
      end
    end
  elseif !circular_variation && extreme_variation  # randomly choose between low_value and high_value
    # TODO:  move these parameters to the configuration file
    # Values of 3/15
    #low_value = 0.1
    #high_value = 0.9
    #range = 0.1
    # values of 3/16
    low_value = 0.25
    high_value = 0.75
    range = 0.02
    for j = 1:num_subpops
      for k = 1:num_attributes
        if rand() < 0.5
          ideal_properties[j].ideal[k] = low_value+(rand()*range-0.5*range)
        else 
          ideal_properties[j].ideal[k] = high_value+(rand()*range-0.5*range)
        end
      end
      #println("j: ",j,"  ideal: ",ideal_properties[j].ideal)
    end
  end
  return ideal_properties
end  

@doc """ horiz_transfer_circular!()
  Transfers variants between subpopulations in a circular fashion (either forward or backward).
  Elements to be transfered are selected by proportional selection.
  Elements to be replaced can be random or selected by reverse proportional selection depending on the flag neg_select.
  subpops is modified by this function (as a side effect)
"""
function horiz_transfer_circular!( N::Int64, m::Int64, num_emmigrants::Int64, subpops::PopList, variant_table::Dict{Int64,variant_type};
     forward::Bool=true, neg_select::Bool=true, emmigrant_select::Bool=true )
  n = Int(floor(N/m))    # size of subpopulations
  #println("horiz_transfer_circular! forward: ",forward)
  emmigrants = PopList()
  for j = 1:m
    if emmigrant_select
      Base.push!( emmigrants, propsel( subpops[j], num_emmigrants, variant_table ) )
    else
      Base.push!( emmigrants, subpops[j][1:num_emmigrants] )   # Neutral
    end
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
  global ideal_properties = ideal_type[ ideal_type(0.1,zeros(Float64,m)) ]
  #global idfe() = 1.0+rand(Distributions.Gamma(1.0,0.1))
  #global cdfe() = -rand(Distributions.Gamma(0.2,0.001))
  global N = 20
  global m = 4
  global mu = 0.2
  global cperr = 0.2
  global ne = 3  # num_emmigrants 
  global ngens = 4
end  
# Sample call to main function
# spatial_simulation( N, m, mu, cperr, ngens, ne, cdfe, idfe, vtbl, ideal_properties )
