#=
Stores the properties of an "innovation" (possibly deleterious or advantageous).
=#
export innovation_type, innovation,  iupdate!, make_extinct!,  make_fixed!
using Base.Test


type innovation_type
  identifier::Int64   # integer for this innovation, not sure if this is needed
  start_gen::Int64    # the generation (time step) when the innovation was generated
  final_gen::Int64    # the generation when the innovation went extinct.  Should be zero while innovation is evolving
  fitness_coefficient::Float64    # The selection coefficient of the innovation (allele)
  sum_counts::Int64   # the current accumulated sum of the counts per generation for infinite sites
  sum_heteroz::Float64  # the current accumulated sum of heterozygosities for infinite sites
  previous_allele_freq::Int64  # Allele frequence of the previous generation
  #history::Vector{Int64}   # history[i] = number of copies in the (i-1)th generation after start_gen
end

@doc """ function innovation() 
  Initializes the properties of the innovation.
"""
function innovation( id::Int64, N::Int64, start_gen::Int64, fitness_coef::Float64=1.0 )
  #println("innovation: id: ",id,"  N: ",N,"  start_gen: ",start_gen)
  initial_heteroz = 1.0 - watterson_homozygosity([N-1,1])
  #return innovation_type( id, start_gen, 0, fitness_coef, 1, initial_heteroz, Int64[1] )
  return innovation_type( id, start_gen, 0, fitness_coef, 1, initial_heteroz, 0  )
end

function iupdate!( innov::innovation_type, N::Int64, generation::Int64, new_allele_freq::Int64 )
  #println("iupdate! site: ",innov.identifier,"   prev: ",innov.previous_allele_freq,"  new: ",new_allele_freq)
  #@test generation == length(innov.history) + innov.start_gen 
  if new_allele_freq > 0
    innov.sum_counts += new_allele_freq
    innov.sum_heteroz += 1.0 - watterson_homozygosity([N-new_allele_freq, new_allele_freq])
  end
  #Base.push!( innov.history, new_allele_freq )
  #println("inn update! index: ",innov.identifier,"  gen: ",generation,"  len(hist): ",length(innov.history), "  start_gen: ",innov.start_gen)
  #println("new_allele_freq: ",new_allele_freq,"  sum_counts: ",innov.sum_counts,"  sum_heteroz: ",innov.sum_heteroz)
  innov.previous_allele_freq = new_allele_freq  # save for the next call to update_selected.
  new_allele_freq
end

# Should be called after, but on the same generation, as update!
function make_extinct!( innov::innovation_type, generation::Int64 )
  #println("make extinct")
  #@test generation == length(innov.history) + innov.start_gen - 1
  innov.final_gen = generation
  #Base.push!( innov.history, 0 )
  #innov.history = Vector{Int64}[]
end

# Should be called after, but on the same generation, as update!
function make_fixed!( innov::innovation_type, generation::Int64 )
  #println("make fixed")
  #@test generation == length(innov.history) + innov.start_gen - 1
  if innov.final_gen == 0   # do not set to be fixed more than once
    innov.final_gen = generation
  end
  #Base.push!( innov.history, allele_freq )
  #innov.history = Vector{Int64}[]
end
