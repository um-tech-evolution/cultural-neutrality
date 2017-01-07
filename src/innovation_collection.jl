#=
This would be a class defintion if Julia was an object-oriented language
Stores a collection of innovations.
Innovations are partitioned into 3 subsets, active, fixed, and extinct.
=#

type innovation_collection
  list::Vector{innovation_type} #    Stores a collection of innovations
  active::IntSet   # indices of actively evolving innovations
  fixed::IntSet   # indices of innovations that have fixed
  extinct::IntSet   # indices of innovations that have gone extinct
end

# Constructor for a new empty innovation collection
function innovation_collection( )  
  innovation_collection( innovation_type[], IntSet(), IntSet(), IntSet() )
end

@doc """ push!()
 Adds an innovation to the collection.
"""
function push!( innov_collection::innovation_collection, innov::innovation_type )
  Base.push!( innov_collection.list, innov )
  if innov.final_gen == 0
    Base.push!( innov_collection.active, length(innov_collection.list) )
  elseif innov.history[end] > 0
    Base.push!( innov_collection.fixed, length(innov_collection.list) )
  else
    Base.push!( innov_collection.extinct, length(innov_collection.list) )
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
  for index in ic.active  # updates sites to the next generation
    #println("innovation: start gen: ",ic.list[index].start_gen,"  history: ",ic.list[index].history)
    new_allele_freq = update_neutral( index, N, ic.list[index].history[end] )
    #println("new_allele_freq: ",new_allele_freq)
    update!(ic,index,g,new_allele_freq)
    if new_allele_freq == 0  # extinction
      make_extinct!(ic,index,g)
    elseif fixation_test( N, new_allele_freq )  # fixation
      make_fixed!(ic,index,g)
    end
  end
end

# Update all active innovations, and make some extinct and fixed
# popcounter is a dictionary that maps alleles to their frequencies in the current generation population
function update_innovations!( ic::innovation_collection, g::Int64, N::Int64, popcounter::Dict{Int64,Int64}, fixation_test::Function=fix_test )
  for index in ic.active  # updates sites to the next generation
    println("update_innovations! index: ",index)
    new_allele_freq = get( popcounter, index, 0 )
    update!(ic,index,g,new_allele_freq)
    if new_allele_freq == 0  # extinction
      make_extinct!(ic,index,g)
    elseif fixation_test( N, new_allele_freq )  # fixation
      make_fixed!(ic,index,g)
    end
  end
end  

# Calls the update! function of innovation on innov_collection.list[index]
function update!( innov_collection::innovation_collection, index::Int64, generation::Int64, new_allele_freq::Int64 )
  #println("i_c update! new_allele_freq: ",new_allele_freq)
  #@test generation == length(innov_collection.list[index].history) + innov_collection.list[index].start_gen
  update!( innov_collection.list[index], generation, new_allele_freq )
end

function make_extinct!( innov_collection::innovation_collection, index::Int64, generation::Int64 )
  make_extinct!( innov_collection.list[index], generation )
  Base.pop!( innov_collection.active,index)
  Base.push!( innov_collection.extinct,index)
end

function make_fixed!( innov_collection::innovation_collection, index::Int64, generation::Int64 )
  make_fixed!( innov_collection.list[index], generation )
  Base.pop!( innov_collection.active,index)
  Base.push!( innov_collection.fixed,index)
end

function average_time_to_extinction( innov_collection::innovation_collection )
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
  if length(innov_collection.fixed) == 0
    println("no fixed innovations")
    return -1.0
  end
  sum = 0
  for i in innov_collection.fixed 
    su += length(innov_collection.list[i].history)
  end
  return Float64(sum-1)/length(innov_collection.fixed)
end 

# Fraction of fixed innovations out of fixed plus extinct innovations
function fixed_fraction( innov_collection::innovation_collection )
  if  length(innov_collection.extinct) == 0 && length(innov_collection.fixed) == 0
    error("no extinct and no fixed innovations")
  end
  return length(innov_collection.fixed)/Float64(length(innov_collection.fixed)+length(innov_collection.extinct))
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
  return ic
end

# TODO:  move to test directory
# Tests some of the functionality of  innovation.jl and innovation_collection.jl
function testall( N::Int64, ngens::Int64, mu::Float64, prob_ext::Float64, prob_fix::Float64 )
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
    done = (g > ngens) && length(ic.active) == 0
  end
  println("number active: ",length(ic.active))
  println("number fixed: ",length(ic.fixed))
  println("number extinct: ",length(ic.extinct))
  println("avg time to fixation: ",average_time_to_fixation(ic))
  println("avg time to extinction: ",average_time_to_extinction(ic))
  println("fixed fraction: ",fixed_fraction(ic))
  ic
end

#=
include("innovation.jl")
include("innovation_collection.jl")
=#
