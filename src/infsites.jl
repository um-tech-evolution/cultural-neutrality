#export infs_result_type, infs_result, site_type, inf_sites, update_neutral, update_selected
export site_type, inf_sites, update_neutral, update_selected
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

type site_type
  id::Int64     # the unique id of this site (i. e., trait or allele)
  count::Int64  # the curent number of copies
  sel_coef::Float64  # selection coefficient where 1 is neutral.   Does not change
end


function push!( innovation::innovation_type, allele_freq::Int64 )
  Base.push!( innovation.history, allele_freq )
end

@doc """ function neutral_inf_sites( )
Do a neutral infinite sites simulation with popsize N, L loci, mutation rate mu per locus,
    ngens generations, and a burn in time of burn_in*N/N_mu+50 generations.
Returns a infs_result_type object, but also prints relevant information.
"""
function inf_sites( N::Int64, N_mu::Float64, ngens::Int64; 
    dfe::Function=dfe_neutral, burn_in::Float64=2.0)
    #ic::innovation_collection=innovation_collection(false) )
  fitness_table = Dict{Int64,Float64}()
  int_burn_in = Int(round(burn_in*N/N_mu+50.0))
  ic = innovation_collection( 1.0 )   # fix_minimum = 1.0 so fixes only if count is N
  g_limit = 1000 
  mu = N_mu/N
  #println("N: ",N)
  #println("N_mu: ",N_mu)
  #println("ngens: ",ngens)
  site_list = 
  poplist= Population[]  # First population has no sites
  new_id = 1
  i = 1
  N_sum = 0
  g = 2
  done = false
  while !done && g < g_limit+ngens+int_burn_in
    #println("generation: ",g)
    update_innovations!( ic, g, N )
    num_mutations = rand(Binomial(N,mu))
    #println("g: ",g,"  num_mutations: ",num_mutations)
    if g < ngens + int_burn_in
      for j in 1:num_mutations
        fit = dfe_fitness( new_id, dfe, fitness_table )  # Set fitness of new_id
        ic_push!(ic,innovation(new_id,g,fit))  
        #println("new innovation: start gen: ",ic.list[end].start_gen,"  history: ",ic.list[end].history)
        new_id += 1
      end
    end
    g += 1
    done = (g > ngens) && length(ic.active) == 0
    N_sum += N_inf_sites(ic)
    #print(" ",N_inf_sites(ic))
  end
  println("Naverage: ", Float64(N_sum)/ngens )
  #print_summary( ic )
  ic
end

@doc """ function update_neutral() 
  Updates site according to the Wright-Fisher model of drift.  
  Since mutation has already happened by the time a new site is established, there is no mutation.
""" 
function update_neutral( site::Int64, N::Int64, old_allele_freq::Int64 )
  p = Float64(old_allele_freq)/N
  new_allele_freq = rand(Binomial(N,p))
  #println("update_neutral: ",site,"  N: ",N,"  new_allele_freq: ",new_allele_freq)
  return new_allele_freq
end 

@doc """ function update_selected() 
  Updates site according to the Wright-Fisher model of drift with selection.  
  Since mutation has already happened by the time a new site is established, there is no mutation.
""" 
function update_selected( site::Int64, N::Int64, old_allele_freq::Int64, select_coef::Float64 )
  p = min(1.0,(select_coef*Float64(old_allele_freq)/N))
  new_allele_freq = rand(Binomial(N,p))
  #println("update_neutral: ",site,"  N: ",N,"  new_allele_freq: ",new_allele_freq)
  return new_allele_freq
end  

# Set some values for testing
#=
n=N=10
ngens=9
N_mu = 1.0
burn_in = 0.0
=#
