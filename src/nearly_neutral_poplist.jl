export fitness, nearly_neutral_poplist, dfe_deleterious,dfe_advantageous, dfe_mixed

using Distributions

@doc """ function fitness( p::Int64, dfe::Function )
Fitness function as saveed in dictionary fitness_table.
"""
function fitness( p::Int64, dfe::Function, nn_select::Int64 )
  global fitness_table = Dict{Int64,Float64}()
  val = get(fitness_table,p,-1.0)
  if val == -1.0   # p is not in fitness table
    val = dfe( p, nn_select )
    fitness_table[p] = val
  end
  val
end 

@doc """ function nearly_neutral_poplist( N::Int64, mu::Float64, ngens::Int64, dfe::Function; burn_in::Float64=1.0,
    uniform_start::Bool=false )
Note:  dfe is "distribution of fitness effects" vector.  
"""
function nearly_neutral_poplist( N::Int64, N_mu::Float64, ngens::Int64, dfe::Function; burn_in::Float64=1.0, 
    uniform_start::Bool=false, nn_select::Int64=1, combine::Bool=true )
  int_burn_in = Int(round(N*burn_in))
  mu = N_mu/N
  if uniform_start  # All allele values start with the same value.  Start with a selective sweep.
    poplist= Population[ Int64[1 for i = 1:N] ]
    fitness(1, dfe, 0 )  # set fitness of 1 to be 1.0
    new_id = 2
  else
    poplist= Population[ collect(1:N) ]
    for i = 1:N
      fitness(i, dfe, 0 )  # set fitness of i to be 1.0
    end
    new_id = N+1
  end
  if combine
    pop_result = Population()
  end
  for g = 2:(ngens+int_burn_in)
    for i in 1:N
      if rand() < mu
        poplist[g-1][i] = new_id
        fitness( new_id, dfe, nn_select )  # Set fitness of new_id
        new_id += 1
      end
    end
    new_pop = propsel( poplist[g-1], fitness, dfe, nn_select )
    push!( poplist, new_pop )
    if combine && g >= int_burn_in+1
      pop_result = vcat( pop_result, new_pop )
    end
  end
  if combine
    return pop_result
  else
    return poplist[int_burn_in+1:end]
  end
end

function dfe_deleterious( x::Int64, nn_select::Int64; alpha::Float64=0.2, theta::Float64=0.5 )
  global dist_deleterious
  if !isdefined(:dist_deleterious)
    dist_deleterious = Distributions.Gamma(alpha,theta)
  end
  if nn_select == 1
    return 1.0-rand(dist_deleterious)
  else
    return 1.0  # flat fitness for no selection
  end
end
  
function dfe_advantageous( N::Int64, nn_select::Int64; alpha::Float64=1.0, theta::Float64=0.01 )
  global dist_advantageous
  if !isdefined(:dist_advantageous)
    dist_advantageous = Distributions.Gamma(alpha,theta)
  end
  if nn_select == 1
    return 1.0+rand(dist_advantageous)
  else
    return 1.0  # flat fitness for no selection
  end
end

function dfe_mixed( N::Int64, nn_select::Int64; adv_probability::Float64=0.2, alpha_disadv::Float64=0.2, alpha_adv::Float64=1.0, theta_disadv::Float64=1.0, theta_adv::Float64=0.01 )
  global dist_deleterious
  if !isdefined(:dist_deleterious)
    dist_deleterious = Distributions.Gamma(alpha_disadv,theta_disadv)
  end
  global dist_advantageous
  if !isdefined(:dist_advantageous)
    dist_advantageous = Distributions.Gamma(alpha_adv,theta_adv)
  end
  if nn_select == 1
    rand() < adv_probability ? 1.0+rand(dist_advantageous) : 1.0-rand(dist_disadvantageous)
  else
    return fill(1.0,N)  # flat fitness for no selection
  end
end
