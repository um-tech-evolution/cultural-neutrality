#=
A simplified version of the infinite sites model where sites evolve independently.
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
    ngens generations, and a burn in time of burn_in*N/N_mu+50 generations.
Returns a sim_result_type object, but also prints relevant information.
"""
function neutral_inf_sites( N::Int64, L::Int64, N_mu::Float64, ngens::Int64; burn_in::Float64=2.0,
    ic::innovation_collection=innovation_collection(false) )
  g_limit = 1000
  int_burn_in = Int(round(burn_in*N/N_mu+50.0))
  mu = N_mu/N
  println("N: ",N)
  println("N_mu: ",N_mu)
  println("ngens: ",ngens)
  if ic.in_use
    println("fix minimum: ",ic.fix_minimum)
  end
  poplist= Population[ collect(1:N) ]
  new_id = 1
  for i = 1:N
    # Note that dfe_advantageous(i), dfe_deleterious(i), dfe_mixed(i) do not depend on i, but dfe_mod(i) does depend on i
    fit = dfe_fitness(i, dfe, fitness_table )  # set fitness of i to be dfe(i).
    new_id += 1
  end
  if combine
    pop_result = Population()
  end
  i = 1
  g = 2
  done = false
  while !done && g < g_limit+ngens+int_burn_in
    #println("generation: ",g)
    update_innovations!( ic, g, N )
    #num_mutations = rand(Poisson(mu*L*N))
    #println("g: ",g,"  num_mutations: ",num_mutations)
    #sim_result.number_mutations += num_mutations
    for j in 1:num_mutations
      push!(ic,innovation(i,g))
      #println("new innovation: start gen: ",ic.list[end].start_gen,"  history: ",ic.list[end].history)
      i += 1
    end
    for i in 1:N
      if rand() < mu
        poplist[g-1][i] = new_id
        fit = dfe_fitness( new_id, dfe, fitness_table )  # Set fitness of new_id
        if g > int_burn_in && g <= ngens+int_burn_in
          #println("id: ",new_id,"  fit: ",fit)
          ic_push!( ic, innovation( new_id, g, fit ) )
        end
        new_id += 1
      end
    end
    g += 1
    done = (g > ngens) && length(ic.active) == 0
  end
  print_summary( ic )
  ic
end

@doc """ function update_site( site::Int64, sim_result::sim_result_type )
  Updates site according to the Wright-Fisher model of drift.  
  Since mutation has already happened by the time a new site is established, there is no mutation.
""" 
function update_neutral( site::Int64, N::Int64, old_allele_freq::Int64 )
  p = Float64(old_allele_freq)/N
  new_allele_freq = rand(Binomial(N,p))
  #println("update_neutral: ",site,"  N: ",N,"  new_allele_freq: ",new_allele_freq)
  return new_allele_freq
end 

function update_selected( site::Int64, N::Int64, old_allel_freq::Int64, select_coef::Float64 )
  p = Float64(N-old_allele_freq)/(N+select_coef*old_allele_freq)
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
