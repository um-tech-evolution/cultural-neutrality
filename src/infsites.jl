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
=#

#using Distributions


@doc """ function neutral_inf_sites( )
Do a neutral infinite sites simulation with popsize N, L loci, mutation rate mu per locus,
    ngens generations, and a burn in time of burn_in*N/N_mu+50 generations.
Returns a infs_result_type object, but also prints relevant information.
"""
function inf_sites( N::Int64, N_mu::Float64, ngens::Int64; 
    dfe::Function=dfe_neutral, burn_in::Float64=2.0)
    #ic::innovation_collection=innovation_collection(N,false) )
  fitness_table = Dict{Int64,Float64}()
  int_burn_in = Int(round(burn_in*N/N_mu+50.0))
  ic = innovation_collection( N, 1.0 )   # fix_minimum = 1.0 so fixes only if count is N
  g_limit = 1000 
  mu = N_mu/N
  #println("N: ",N)
  #println("N_mu: ",N_mu)
  #println("ngens: ",ngens)
  #site_liss_table = Dict{Int64,Float64}()
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
        ic_push!(ic,innovation(new_id,N,g,fit))  
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
  #(count_adv_fixed, count_del_fixed) = count_adv_del_fixed( ic )
  #ic.count_fixed_adv = count_adv_fixed
  #ic.count_fixed_del = count_del_fixed
  (ic.count_fixed_adv, ic.count_fixed_del) = count_adv_del_fixed( ic )
  println("count_fixed_del: ",ic.count_fixed_del)
  println("count_fixed_adv: ",ic.count_fixed_adv)
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
  This function is called by update_innovations!() in innovation_collection.jl.

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
