#=
Neutral infinite sites simulation where sites evolve independently.
 N = popsize
 n = sample size
 L = sequence length
 mu = per site mutation rate
 ngens = number generations
 burn_in
TODO:  keep track of the locus, varying from 1 to L, for each mutation.  
     This is different from the site used to simulate what happens to the mutation.
     The infinite sites model assumes that each mutation is at a new site.  
     Thus, a mutation could happen at a locus where there is already an active site----
      this would be a case where the model assumptions were violated.  The program would
      create a new site to keep track of the mutation independently of the previous site 
      corresponding to this locus.
=#

using Distributions

#=
type innovation_type
  selection_coefficient::Float64    # The selection coefficient of the innovation (allele)
  start_gen::Int64    # the generation (time step) when the innovation was generated
  final_gen::Int64    # the generation when the innovation went extinct
  history::Vector{Int64}   # history[i] = number of copies in the (i-1)th generation after start_gen
end
=#

type sim_result_type
  N::Int64
  L::Int64
  mu::Float64
  ngens::Int64
  number_mutations::Int64
  number_extinctions::Int64
  number_fixations::Int64
  total_extinction_time::Int64
  total_fixation_time::Int64
  sites::Vector{innovation_type}   # includes the allele frequency as it varies over time
end

function push!( innovation::innovation_type, allele_freq::Int64 )
  Base.push!( innovation.history, allele_freq )
end

@doc """ function neutral_inf_sites( N::Int64, L::Int64, mu::Float64, ngens::Int64, burn_in::Float64 )
Do a neutral infinite sites simulation with popsize N, L loci, mutation rate mu per locus,
    ngens generations, and a burn in time of burn_in*N generations.
Returns a sim_result_type object, but also prints relevant information.
"""
function neutral_inf_sites( N::Int64, L::Int64, mu::Float64, ngens::Int64, burn_in::Float64 )
  int_burn_in = Int(round(N*burn_in))
  #global innovation_list = innovation_type[]
  ic = innovation_collection()
  #sim_result = sim_result_type( N, L, mu, ngens, 0,0,0,0,0, locus_type[] )
  i = 1
  g = 1
  g_limit = 1000
  done = false
  while !done && g < g_limit
    println("generation: ",g)
    update_innovations!( ic, g, N )
    #=
    for index in ic.active  # updates sites to the next generation
      println("innovation: start gen: ",ic.list[index].start_gen,"  history: ",ic.list[index].history)
      new_allele_freq = update_neutral( index, N, ic.list[index].history[end] )
      println("new_allele_freq: ",new_allele_freq)
      update!(ic,index,g,new_allele_freq)
      if new_allele_freq == 0  # extinction
        make_extinct!(ic,index,g)
      elseif new_allele_freq == N  # fixation
        make_fixed!(ic,index,g)
      end
    end
    =#
    num_mutations = rand(Poisson(mu*L*N))
    println("g: ",g,"  num_mutations: ",num_mutations)
    #sim_result.number_mutations += num_mutations
    for j in 1:num_mutations
      push!(ic,innovation(i,g))
      #println("new innovation: start gen: ",ic.list[end].start_gen,"  history: ",ic.list[end].history)
      i += 1
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

@doc """ function update_site( site::Int64, sim_result::sim_result_type )
  Updates site according to the Wright-Fisher model of drift.  
  Mutation is handled elsewhere.
""" 
function update_neutral( site::Int64, N::Int64, old_allele_freq::Int64 )
  p = Float64(old_allele_freq)/N
  new_allele_freq = rand(Binomial(N,p))
  #println("update_neutral: ",site,"  N: ",N,"  new_allele_freq: ",new_allele_freq)
  return new_allele_freq
end 
  
# Set some values for testing
n=N=10
ngens=9
mu = 1.0/N
burn_in = 0.0
L=2
  
