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

function time_series( n::Int64, N::Int64, mu::Float64, ngens::Int64, Btbl::Array{Float64,2}, burn_in::Float64=2.0 )
  pop_list = neutral_poplist(N, mu, ngens, burn_in=burn_in )
  p_counts_list = map(x->pop_counts64(sample_population(x,n)),pop_list)
  cfg_list = Any[]
  w_heteroz_list = Float64[]
  s_prob_list = Float64[]
  e_prob_list = Float64[]
  ewens_theta_list = Float64[]
  for cfg in p_counts_list
    w_heteroz = 1.0 - watterson_homozygosity(cfg)
    #s_enum = slatkin_enum( cfg )  # not working correctly
    ewens_m = ewens_montecarlo(100000,cfg)
	  ucfg =  ConfigInt[ x for x in cfg ]
	  #println("ucfg: ",ucfg)
	  s_prob = slatkin_exact( ucfg, Btbl )
	  #println("s_prob: ",s_prob)
    push!( cfg_list, cfg )
    push!( w_heteroz_list, w_heteroz )
    push!( s_prob_list, s_prob )
    push!( e_prob_list, ewens_m.probability )
    push!( ewens_theta_list, ewens_m.theta_estimate )
  end
  DataFrame(
    cfg = cfg_list,
    w_heteroz = w_heteroz_list,
    s_prob = s_prob_list,
    e_prob = s_prob_list,
    ewens_theta = ewens_theta_list,
  )
end

function test_independence( ts_df::DataFrame, sym::Symbol, sig_level::Float64, incr::Int64=1 )
	contig_table = zeros(Float64, 3, 3)
  expected_table = zeros(Float64, 3, 3)
	for i = 1:(size(ts_df)[1]-incr)
		L1 =  ts_df[sym][i] < sig_level 
    L2 =  ts_df[sym][i+incr] < sig_level
		G1 =  ts_df[sym][i] > 1.0 - sig_level 
		G2 =  ts_df[sym][i+incr] > 1.0 - sig_level 
		contig_table[1,1] += L1 && L2 ? 1 : 0
		contig_table[1,2] += L1 && !L2 && !G2 ? 1 : 0
		contig_table[1,3] += L1 && G2 ? 1 : 0
		contig_table[2,1] += !L1 && !G1 && L2 ? 1 : 0
		contig_table[2,2] += !L1 && !G1 && !L2 && !G2 ? 1 : 0
		contig_table[2,3] += !L1 && !G1 && G2 ? 1 : 0
		contig_table[3,1] += G1 && L2 ? 1 : 0
		contig_table[3,2] += G1 && !L2 && !G2 ? 1 : 0
		contig_table[3,3] += G1 && G2 ? 1 : 0
	end
  sum_contig = sum(contig_table)
  contig_table /= sum_contig 
  sum1 = [ sum(contig_table[:,i]) for i = 1:3 ]
  sum2 = [ sum(contig_table[i,:]) for i = 1:3 ]
  for i = 1:3
    for j = 1:3
      expected_table[i,j] = sum1[i]*sum2[j]
    end
  end
	contig_table, expected_table
end	

function slat_test(ts_df::DataFrame, sym::Symbol, sig_level::Float64 )
  count_less = 0.0
  count_more = 0.0
  count_mid = 0.0
  nreps = size(ts_df)[1]
	for i = 1:nreps
		S1 =  ts_df[sym][i] < sig_level 
    S3 =  ts_df[sym][i] > 1.0 - sig_level
    S2 =  !S1 && !S3
	  count_less += S1 ? 1.0 : 0.0	
	  count_mid += S2 ? 1.0 : 0.0	
	  count_more += S3 ? 1.0 : 0.0	
  end  
  ( count_less/nreps, count_mid/nreps, count_more/nreps )
end
		
    
function acor(lst, lags::Vector{Int64}=collect(1:10))
	StatsBase.autocor(Float64[x for x in lst],lags)
end

function acor_ts( n::Int64, N::Int64, mu::Float64, ngens::Int64, burn_in::Float64=2.0 )
	ts = time_series( n, N, mu, ngens, burn_in )
  DataFrame(
		w_het_acor = acor( ts[:w_heteroz] ),
		s_prob_acor = acor( ts[:s_prob] )
	)
end

function mult_acor_ts( nts::Int64, n::Int64, N::Int64, mu::Float64, ngens::Int64, burn_in::Float64=2.0 )
  mts = pmap( x->acor_ts(n,N,mu,ngens,burn_in), fill(1,nts))
  DataFrame(
		ave_whet_acor = [ mean( [ mts[i][:w_het_acor][j] for i in collect(1:nts) ]) for j in collect(1:10) ],
		ave_whet_sprob = [ mean( [ mts[i][:s_prob_acor][j] for i in collect(1:nts) ]) for j in collect(1:10) ]
  )
end
