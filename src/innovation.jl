#=
Code for an innovation class that stores an innovation (mutation).
=#
using Base.Test

type innovation_type
  identifier::Int64   # integer for this innovation, not sure if this is needed
  start_gen::Int64    # the generation (time step) when the innovation was generated
  final_gen::Int64    # the generation when the innovation went extinct.  Should be zero while innovation is evolving
  selection_coefficient::Float64    # The selection coefficient of the innovation (allele)
  history::Vector{Int64}   # history[i] = number of copies in the (i-1)th generation after start_gen
end

@doc """ function innovation() 
  Constructor for an innovation object.
"""
function innovation( id::Int64, start_gen::Int64, selection_coef::Float64=1.0 )
  return innovation_type( id, start_gen, 0, selection_coef, Int64[1] )
end

function update!( innov::innovation_type, generation::Int64, new_allele_freq::Int64 )
  println("inn update! index: ",innov.identifier,"  gen: ",generation,"  len(hist): ",length(innov.history),
    "  start_gen: ",innov.start_gen)
  @test generation == length(innov.history) + innov.start_gen 
  Base.push!( innov.history, new_allele_freq )
  new_allele_freq
end

# Should be called after, but on the same generation, as update!
function make_extinct!( innov::innovation_type, generation::Int64 )
  @test generation == length(innov.history) + innov.start_gen - 1
  innov.final_gen = generation
  #Base.push!( innov.history, 0 )
end

# Should be called after, but on the same generation, as update!
function make_fixed!( innov::innovation_type, generation::Int64 )
  @test generation == length(innov.history) + innov.start_gen - 1
  if innov.final_gen == 0   # do not set to be fixed more than once
    innov.final_gen = generation
  end
  #Base.push!( innov.history, allele_freq )
end
