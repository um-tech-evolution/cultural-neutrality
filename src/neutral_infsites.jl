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

type locus_type
  k::Int64   # current number of copies of allele at this site
  start_gen::Int64    # the generation (time step) when the allele was generated
end

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
  sites::Vector{locus_type}   # includes the allele frequency as it varies over time
end
@doc """ function neutral_inf_sites( N::Int64, L::Int64, mu::Float64, ngens::Int64, burn_in::Float64 )
Do a neutral infinite sites simulation with popsize N, L loci, mutation rate mu per locus,
    ngens generations, and a burn in time of burn_in*N generations.
Returns a sim_result_type object, but also prints relevant information.
"""
function neutral_inf_sites( N::Int64, L::Int64, mu::Float64, ngens::Int64, burn_in::Float64 )
  int_burn_in = Int(round(N*burn_in))
  num_sites = 0  # number of sites in use
  next_avail_site = 1
  active_sites = IntSet()  # Sites that are in use (polymorphic)
  avail_sites = IntSet() # Already allocated sites that are not in use
  
  sim_result = sim_result_type( N, L, mu, ngens, 0,0,0,0,0, locus_type[] )
  for g = 1:int_burn_in+ngens
    for key in active_sites  # updates sites to the next generation
      new_k = update_site( key, sim_result )
      if new_k == 0  # extinction
        push!(avail_sites,key)
        pop!(active_sites,key)
        extinction_time = g - sim_result.sites[key].start_gen
        #println("extinction at site: ",key," after ", extinction_time," generations")
        sim_result.number_extinctions += 1
        sim_result.total_extinction_time += extinction_time
      elseif new_k == N  # fixation
        push!(avail_sites,key)
        pop!(active_sites,key)
        fixation_time = g - sim_result.sites[key].start_gen
        #println("fixation at site: ",key," after ", fixation_time," generations")
        sim_result.number_fixations += 1
        sim_result.total_fixation_time += fixation_time
      end
    end
    num_mutations = rand(Poisson(mu*L*N))
    #println("g: ",g,"  num_mutations: ",num_mutations)
    sim_result.number_mutations += num_mutations
    for i in 1:num_mutations
      if length(avail_sites) > 0
        new_site = pop!(avail_sites)
        #println("new site: ",new_site)
        push!(active_sites,new_site) 
        sim_result.sites[new_site].k = 1 # Initialize new site with a allele count of 1
        sim_result.sites[new_site].start_gen = g # Initialize new site with a allele count of 1
      else
        new_site = length(sim_result.sites)+1
        #println("new site: ",new_site)
        push!(active_sites,new_site)
        push!(sim_result.sites,locus_type(1,g))   # Initialize new site with a trait count of 1
      end
      #=
      println("sites: ")
      for i = 1:length(sim_result.sites)
        println("sites[",i,"] = ",sim_result.sites[i])
      end
      =#
    end
  end
  println("number extinctions: ",sim_result.number_extinctions)
  println("average extinction time: ",
      Float64(sim_result.total_extinction_time)/sim_result.number_extinctions)
  println("number fixations: ",sim_result.number_fixations)
  println("average fixation time: ",
      Float64(sim_result.total_fixation_time)/sim_result.number_fixations)
  println("number sites: ",length(sim_result.sites))
  println("number active sites: ",length(active_sites))
  println("number avail sites: ",length(avail_sites))
  return sim_result
end

@doc """ function update_site( site::Int64, sim_result::sim_result_type )
  Updates site according to the Wright-Fisher model of drift.  
  Mutation is handled elsewhere.
""" 
function update_site( site::Int64, sim_result::sim_result_type )
  k = sim_result.sites[site].k
  p = Float64(k)/sim_result.N
  new_k = rand(Binomial(N,p))
  sim_result.sites[site].k = new_k
  #println("update_site: ",site,"  new_k: ",new_k)
  return new_k
end 
  
# Set some values for testing
n=N=10
ngens=16
mu = 3.0/N
burn_in = 0.0
L=3
  
