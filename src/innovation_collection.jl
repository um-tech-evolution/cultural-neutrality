#=
This would be a class defintion if Julia was an object-oriented language
Stores a collection of innovations.
Innovations are partitioned into 3 subsets, active, fixed, and extinct.
=#

export innovation_collection, ic_push!, update_innovations!, ic_update!, make_extinct!, 
    make_fixed!, print_summary, average_time_to_extinction, average_time_to_fixation, 
    fixed_fraction, average_fitness_fixed, average_fitness_extinct, average_fitness_all, fix_test

type innovation_collection
  list::Dict{Int64,innovation_type} #    Stores a collection of innovations
  active::IntSet   # indices of actively evolving innovations
  fixed::IntSet   # indices of innovations that have fixed
  extinct::IntSet   # indices of innovations that have gone extinct
  fix_minimum::Float64  # minimum fraction of popsize for fixation
  in_use::Bool      # If false, not used
end

# Constructor for a new empty innovation collection
function innovation_collection(in_use::Bool=true )  
  innovation_collection( Dict{Int64,innovation_type}(), IntSet(), IntSet(), IntSet(), 1.0, in_use )
end

# Constructor for a new empty innovation collection with a value for fix_minimum
function innovation_collection( fix_min::Float64, in_use::Bool=true )  
  innovation_collection( Dict{Int64,innovation_type}(), IntSet(), IntSet(), IntSet(), fix_min, in_use )
end

@doc """ ic_push!()
 Adds an innovation to the collection.
"""
function ic_push!( innov_collection::innovation_collection, innov::innovation_type )
  if !innov_collection.in_use 
    return
  end
  #println("new innovation: id: ",innov.identifier,"  sg: ",innov.start_gen,"  fg: ",innov.final_gen,"  len(hist): ",length(innov.history))
  innov_collection.list[innov.identifier] = innov 
  if innov.final_gen == 0
    Base.push!( innov_collection.active, innov.identifier )
  elseif innov.history[end] > 0
    Base.push!( innov_collection.fixed, innov.identifier )
  else
    Base.push!( innov_collection.extinct, innov.identifier )
  end
  innov_collection
end

@doc """ fix_test()
The most basic test for fixation, namely that the allele frequency is the population size.
"""
function fix_test( N::Int64, new_allele_freq::Int64 )
  return N == new_allele_freq
end

# Update all active innovations, and make some extinct and fixed
function update_innovations!( ic::innovation_collection, g::Int64, N::Int64, fixation_test::Function=fix_test )
  if !ic.in_use 
    return
  end
  for index in ic.active  # updates sites to the next generation
    #println("innovation: start gen: ",ic.list[index].start_gen,"  history: ",ic.list[index].history)
    new_allele_freq = update_neutral( index, N, ic.list[index].history[end] )
    #println("new_allele_freq: ",new_allele_freq)
    ic_update!(ic,index,g,new_allele_freq)
    if new_allele_freq == 0  # extinction
      make_extinct!(ic,index,g)
    elseif new_allele_freq >= Int(ceil(ic.fix_minimum*N))  # fixation
      make_fixed!(ic,index,g)
    end
  end
end

# Update all active innovations, and make some extinct and fixed
# popcounter is a dictionary that maps alleles to their frequencies in the current generation population
function update_innovations!( ic::innovation_collection, g::Int64, N::Int64, popcounter::Dict{Int64,Int64}, fixation_test::Function=fix_test )
  if !ic.in_use 
    return
  end
  #println("update_innovations! active: ",ic.active)
  for index in ic.active  # updates sites to the next generation
    #println("update_innovations! index: ",index)
    new_allele_freq = get( popcounter, index, 0 )
    ic_update!(ic,index,g,new_allele_freq)
    if new_allele_freq == 0  # extinction
      make_extinct!(ic,index,g)
    elseif new_allele_freq >= Int(ceil(ic.fix_minimum*N) )  # fixation
      make_fixed!(ic,index,g)
    end
  end
end  

# Calls the iupdate! function of innovation on innov_collection.list[index]
function ic_update!( innov_collection::innovation_collection, index::Int64, generation::Int64, new_allele_freq::Int64 )
  if !innov_collection.in_use 
    return
  end
  #println("i_c update! new_allele_freq: ",new_allele_freq)
  #@test generation == length(innov_collection.list[index].history) + innov_collection.list[index].start_gen
  iupdate!( innov_collection.list[index], generation, new_allele_freq )
end

function make_extinct!( innov_collection::innovation_collection, index::Int64, generation::Int64 )
  if !innov_collection.in_use 
    return
  end
  make_extinct!( innov_collection.list[index], generation )
  Base.pop!( innov_collection.active,index)
  Base.push!( innov_collection.extinct,index)
end

function make_fixed!( innov_collection::innovation_collection, index::Int64, generation::Int64 )
  if !innov_collection.in_use 
    return
  end
  make_fixed!( innov_collection.list[index], generation )
  Base.pop!( innov_collection.active,index)
  Base.push!( innov_collection.fixed,index)
end

# TODO:  Move to conformist_poplist.jl or delete
function compute_turnovers( pop1::Population, pop2::Population, N_mu::Float64, Ylist::Vector{Int64},
    Zsum_list::Vector{Int64}, count_list::Vector{Int64} )
  i = 1
  for y in Ylist
    if Float64(y) < 5.0*N_mu
      Zsum_list[i] += turnover( pop1, pop2, y )
      count_list[i] += 1
    end
    i+= 1
  end
end

function print_summary( ic::innovation_collection; print_lists::Bool=false )
  if !ic.in_use 
    return
  end
  if print_lists
    println("active list: ",ic.active)
    println("fixed list: ",ic.fixed)
    println("extinct list: ",ic.extinct)
  end
  println("number active: ",length(ic.active))
  println("number fixed: ",length(ic.fixed))
  println("number extinct: ",length(ic.extinct))
  println("fixed fraction: ",fixed_fraction(ic))
  println("avg time to fixation: ",average_time_to_fixation(ic))
  println("avg time to extinction: ",average_time_to_extinction(ic))
  println("avg fitness fixed: ",average_fitness_fixed(ic))
  println("avg fitness extinct: ",average_fitness_extinct(ic))
  println("avg fitness all: ",average_fitness_all(ic))
end

function average_time_to_extinction( innov_collection::innovation_collection )
  if !innov_collection.in_use 
    return 0.0
  end
  if length(innov_collection.extinct) == 0
    println("no extinct innovations")
    return -1.0
  end
  sum = 0
  for i in innov_collection.extinct 
    sum += length(innov_collection.list[i].history)
  end
  return Float64(sum-1)/length(innov_collection.extinct)
end 

function average_time_to_fixation( innov_collection::innovation_collection )
  if !innov_collection.in_use 
    return 0.0
  end
  if length(innov_collection.fixed) == 0
    println("no fixed innovations")
    return -1.0
  end
  sum = 0
  for i in innov_collection.fixed 
    sum += length(innov_collection.list[i].history)
  end
  return Float64(sum-1)/length(innov_collection.fixed)
end 

# Fraction of fixed innovations out of fixed plus extinct innovations
function fixed_fraction( innov_collection::innovation_collection )
  if !innov_collection.in_use 
    return 0.0
  end
  if  length(innov_collection.extinct) == 0 && length(innov_collection.fixed) == 0
    error("no extinct and no fixed innovations")
  end
  return length(innov_collection.fixed)/Float64(length(innov_collection.fixed)+length(innov_collection.extinct))
end

function average_fitness_fixed( ic::innovation_collection )
  if !ic.in_use 
    return 0.0
  end
  if  length(ic.fixed) == 0
    return -1.0
  end
  sum = 0.0
  for i in ic.fixed 
    sum += ic.list[i].selection_coefficient
  end
  return sum/length(ic.fixed)
end

function average_fitness_extinct( ic::innovation_collection )
  if !ic.in_use 
    return 0.0
  end
  if  length(ic.extinct) == 0
    return -1.0
  end
  sum = 0.0
  for i in ic.extinct 
    sum += ic.list[i].selection_coefficient
  end
  return sum/length(ic.extinct)
end

function average_fitness_all( ic::innovation_collection )
  if !ic.in_use 
    return 0.0
  end
  if  length(ic.extinct) == 0
    return -1.0
  end
  sum = 0.0
  for i in ic.extinct 
    sum += ic.list[i].selection_coefficient
  end
  for i in ic.fixed 
    sum += ic.list[i].selection_coefficient
  end
  return sum/(length(ic.extinct) + length(ic.fixed))
end

