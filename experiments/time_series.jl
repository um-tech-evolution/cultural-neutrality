using DataFrames
#=  Generate a time series of 4-tuples: 
  1)  Watterson Heterozygosity
  2)  Slatkin test probability
  3)  Watterson test probability
  4)  Slatkin theta estimate
from successive generations of the infinite alleles model after a period of burn-in.
n = sample size
N = populaiton size
mu = mutation (innovation) rate
ngens = number of generations after burn-in
burn_in = multiple of N for burn_in
=#

function time_series( n::Int64, N::Int64, mu::Float64, ngens::Int64, burn_in::Float64 )
  pop_list = neutral_poplist(N, mu, ngens, burn_in=burn_in )
  p_counts_list = map(x->pop_counts64(sample_population(x,n)),pop_list)
  cfg_list = Any[]
  w_heteroz_list = Float64[]
  ewens_m_list = Float64[]
  ewens_theta_list = Float64[]
  for cfg in p_counts_list
    w_heteroz = 1.0 - watterson_homozygosity(cfg)
    #s_enum = slatkin_enum( cfg )  # not working correctly
    ewens_m = ewens_montecarlo(100000,cfg)
    push!( cfg_list, cfg )
    push!( w_heteroz_list, w_heteroz )
    push!( ewens_m_list, ewens_m.probability )
    push!( ewens_theta_list, ewens_m.theta_estimate )
  end
  DataFrame(
    cfg = cfg_list,
    w_heteroz = w_heteroz_list,
    ewens_theta = ewens_theta_list,
    ewens_m = ewens_m_list
  )
end
    
