#=
This would be a class defintion if Julia was an object-oriented language
Stores a collection of innovations.
Innovations are partitioned into 3 subsets, active, fixed, and extinct.
=#

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

@doc """ push!()
 Adds an innovation to the collection.
"""
function push!( innov_collection::innovation_collection, innov::innovation_type )
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
    update!(ic,index,g,new_allele_freq)
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
    update!(ic,index,g,new_allele_freq)
    if new_allele_freq == 0  # extinction
      make_extinct!(ic,index,g)
    elseif new_allele_freq >= Int(ceil(ic.fix_minimum*N) )  # fixation
      make_fixed!(ic,index,g)
    end
  end
end  

# Calls the update! function of innovation on innov_collection.list[index]
function update!( innov_collection::innovation_collection, index::Int64, generation::Int64, new_allele_freq::Int64 )
  if !innov_collection.in_use 
    return
  end
  #println("i_c update! new_allele_freq: ",new_allele_freq)
  #@test generation == length(innov_collection.list[index].history) + innov_collection.list[index].start_gen
  update!( innov_collection.list[index], generation, new_allele_freq )
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


# TODO:  move to test directory
# A "sanity test" for innovation.jl and innovation_collection.jl
function test()
  ic = innovation_collection()
  # Create 2 innovations on generation 1
  push!(ic,innovation(1,1))
  push!(ic,innovation(2,1))
  # Updates for generation 2
  update!(ic,1,2,4)
  update!(ic,2,2,5)
  # Create 1 innovation on generation 2
  push!(ic,innovation(3,2))
  # Make 1 extinct and 2 fixed on generation 3
  update!(ic,1,3,0)
  make_extinct!(ic,1,3)
  update!(ic,2,3,8)
  make_fixed!(ic,2,3)
  # Update 3 on generation 3
  update!(ic,3,3,4)
  # Update 3 on generation 4
  update!(ic,3,4,8)
  make_fixed!(ic,3,4)
  print_summary( ic, print_lists=true )
  return ic
end

# TODO:  move to test directory
# Tests some of the functionality of  innovation.jl and innovation_collection.jl
function testall( N::Int64, ngens::Int64, mu::Float64; prob_ext::Float64=0.3, prob_fix::Float64=0.2 )
  ic = innovation_collection()
  i = 1
  g = 1
  g_limit = 1000
  done = false
  while !done && g < g_limit
    println("generation: ",g)
    for index in ic.active
      r = rand()
      ac = rand(1:N)
      println("updating ",index,"  id:",ic.list[index].identifier)
      update!(ic,index,g,ac)
      if r < prob_ext
        println("extincting ",index,"  id:",ic.list[index].identifier)
        make_extinct!(ic,index,g)
      elseif r < prob_ext+prob_fix
        println("fixing ",index,"  id:",ic.list[index].identifier)
        make_fixed!(ic,index,g)
      end
    end
    if g <= ngens
      for j = 1:N
        r = rand()
        if r < mu
          println("generating innovation ",i)
          push!(ic,innovation(i,g))
          i += 1
        end
      end
    end
    g += 1
    if g == ngens
      println("active: ",ic.active)
    end
    done = (g > ngens) && length(ic.active) == 0
  end
  print_summary( ic, print_lists=true )
  ic
end

function testall()
  testall(N,ngens,mu)
end

#=
include("innovation.jl")
include("innovation_collection.jl")
=#
