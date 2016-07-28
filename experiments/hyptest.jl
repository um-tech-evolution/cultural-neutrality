#=
Analyzes the CSV file produced by running simulation.jl and run.jl.
This CSV file must include the case where cpower == 0 or acer_C == 0 
  which is the unbiased tranmission case.
If the neutral infinite alleles model holds, then the expected homozygosity is 1.0/(1.0+theta)
    where theta = 2*N*mu.
Applies the Watterson and Slatkin homozygosity tests to the null hypothesis that  the homozygosity 
    predicted by the Watterson and Slatkin tests is equal to  1.0/(1.0+theta) where theta = 2*N*mu.
    Note that N and mu are known parameters of the simulation.  
theta == 0.0 against the alternatives:
   theta < 0.0    for conformist transmission
   theta > 0.0    for anti-conformist transmission
Also applies the Slatkin "exact" test (more precicely, the Monte Carlo extension) to the
same null and alternative hypotheses.

In the resulting table,  wq_lo is the wq_lo quantile and wq_hi is the wq_hi quantile.
For the null hypothesis  cpower==0, wq_lo is critcal theta value for the theta<0 alternative hypothesis,
  and wp_hi is the critical theta value for the theta>0 alternative hypothesis.
A type I error is the incorrect rejection of a true null hypothesis, and a type II error is the failure
  to reject a false null hypothesis.
We assume a lo_val type I error probability, and compute the type II error probability using the wq_lo 
  (for conformist) quantile or the wq_hi (for anti-conformist) quantile of the null hypothesis distribution.

The program can be run either as a command-line program, or "included" into Julia.  When run from
the command line, the the "simname" filename parameter must be the command-line argument.
Suggested command-line usage:
 julia -p 8 -L simulation.jl run.jl configs/example1.jl
 julia hyptest.jl configs/example1.jl
Suggested non-command-line usage (assuming configs/example.csv has been computeed by simulation.jl):
 julia
 > include("hyptest.jl")
 > hyptest("configs/example1")
Writes results as a data frame to stdout, and to a CSV file ""$(simname)_hypdf.csv".
=#

using DataFrames
include("dataframe_io.jl")

#TODO:  consider consolidating the "count_more" and "count_less" functions

const lo_val = 0.05         # Tests are at the 5% confidence level
const hi_val = 1.0-lo_val

@doc """ function(count_more_w_homo( df, cvalue::Float64, c_symbol::Symbol, N::Int64, N_mu::Float64, wq_lo::Float64 )

Counts the number of alternative hypothesis outcomes (assuming N, N_mu, cpower as specified)
  which have a w_homo value greater than wq_lo where wq_lo is the lo_val critical value for
  the null hypothesis.
"""
function count_more_w_homo( df, cvalue::Float64, c_symbol::Symbol, N::Int64, N_mu::Float64, wq_lo::Float64 )
  countnz( (df[:N].==N)&(df[:N_mu].==N_mu)&(df[c_symbol].==cvalue)&(df[:w_homoz].>wq_lo))
end

@doc """ function count_less_w_homo( df, cvalue::Float64, c_symbol::Symbol, N::Int64, N_mu::Float64, wq_hi::Float64 )

Counts the number of alternative hypothesis outcomes (assuming N, N_mu, cpower as specified)
  which have a w_homo value less than wq_hi where wq_hi is the lo_val critical value for
  the null hypothesis.
"""
function count_less_w_homo( df, cvalue::Float64, c_symbol::Symbol, N::Int64, N_mu::Float64, wq_hi::Float64 )
  countnz( (df[:N].==N)&(df[:N_mu].==N_mu)&(df[c_symbol].==cvalue)&(df[:w_homoz].<wq_hi))
end

@doc """ function(count_more_s_homo( df, cvalue::Float64, c_symbol::Symbol, N::Int64, N_mu::Float64, sq_lo::Float64 )

Counts the number of alternative hypothesis outcomes (assuming N, N_mu, cpower as specified)
  which have a s_homo value greater than sq_lo where sq_lo is the lo_val critical value for
  the null hypothesis.
"""
function count_more_s_homo( df, cvalue::Float64, c_symbol::Symbol, N::Int64, N_mu::Float64, sq_lo::Float64 )
  countnz( (df[:N].==N)&(df[:N_mu].==N_mu)&(df[c_symbol].==cvalue)&(df[:s_homoz].>sq_lo))
end

@doc """ function count_less_s_homo( df, cvalue::Float64, c_symbol::Symbol, N::Int64, N_mu::Float64, sq_hi::Float64 )

Counts the number of alternative hypothesis outcomes (assuming N, N_mu, cpower as specified)
  which have a s_homo value less than sq_hi where sq_hi is the lo_val critical value for
  the null hypothesis.
"""
function count_less_s_homo( df, cvalue::Float64, c_symbol::Symbol, N::Int64, N_mu::Float64, sq_hi::Float64 )
  countnz( (df[:N].==N)&(df[:N_mu].==N_mu)&(df[c_symbol].==cvalue)&(df[:s_homoz].<sq_hi))
end

@doc """ function s_prob_accept( df, cvalue::Float64, c_symbol::Symbol, N::Int64, N_mu::Float64 )
  Given a cvalue, a value of N, and a value of N_mu, count the number of times that the Slatkin
    test is accepted at the two-sided lo_val level.
"""
function count_s_prob_accept( df, cvalue::Float64, c_symbol::Symbol, N::Int64, N_mu::Float64 )
  countnz( (df[:N].==N)&(df[:N_mu].==N_mu)&(df[c_symbol].==cvalue)&(df[:s_prob].>(lo_val/2.0))&(df[:s_prob].<(1.0-lo_val/2.0)))
end

@doc """ function hyptest_results( df )

Returns a data frame that includes the lo_val and 1.0-lo_val quantiles of the null hypothesis distribution
  and the alternative hypothesis distributions.
Includes zero arrays for the Watterson and Slatkin homozygosity type 2 error and the Slatkin probability type 2 error.
  These will be computered later.
"""
function hyptest_results( df::DataFrame, c_symbol::Symbol )
  result_df = by(df, [ c_symbol, :N_mu, :N ] ) do d
    DataFrame(
      #theta = 
      homozyg = 1.0/(1.0 + 2.0*d[:N_mu][1]),
      mean_w_homo=mean(d[:w_homoz]),
      wq_lo = quantile(d[:w_homoz],lo_val),
      wq_hi = quantile(d[:w_homoz],hi_val),
      w_typ2err = 0.0,
      mean_s_homo=mean(d[:s_homoz]),
      sq_lo = quantile(d[:s_homoz],lo_val),
      sq_hi = quantile(d[:s_homoz],hi_val),
      s_typ2err = 0.0,
      mean_s_prob=mean(d[:s_prob]),
      sp_lo = quantile(d[:s_prob],lo_val),
      sp_hi = quantile(d[:s_prob],hi_val),
      sp_t2err = 0.0,
    )
  end
  #result_df[ df[c_symbol].==0.0 ]
  result_df
end

@doc """ function hyptest( simname )

Adds columns "w_typ2err" (Watterson) and "s_typ2err" (Slatkin)  to the dataframe computed by function 
  hyptest_results().
For the alternative hypothesis rows, this is the type II error corresponding to a type I error
  probability of lo_val1.
"""
function hyptest( simname )
	include("$(simname).jl")
	M = length(N_list)*length(N_mu_list)  # number of cases of N times number of cases of N_mu
	df = readtable("$(simname).csv", makefactors=true, allowcomments=true)
  if findfirst(names(df),:cpower) > 0   # power model of conformity
    c_symbol = :cpower
  elseif findfirst(names(df),:acer_C) > 0  # Acerbi model of conformity
    c_symbol = :acer_C
  else
    error( "file $(simname).csv must contain a column named either :cpower or :acer_C ")
  end
	hypdf = hyptest_results( df, c_symbol )
  sort!(hypdf,cols=[order(c_symbol,by=abs)])
	check_zero_first( M, hypdf, c_symbol )
	w_typ2err = zeros(Float64,size(hypdf,1))
	for i = 1:size(hypdf,1)
    if hypdf[c_symbol][i] > 0.0  # conformity 
      hypdf[:w_typ2err][i] = count_less_w_homo( df, hypdf[c_symbol][i], c_symbol, hypdf[:N][i],hypdf[:N_mu][i],hypdf[:wq_hi][(i-1)%M+1] )/Float64(T)
    elseif hypdf[c_symbol][i] < 0.0  # anti-conformity
      hypdf[:w_typ2err][i] = count_more_w_homo( df, hypdf[c_symbol][i], c_symbol, hypdf[:N][i],hypdf[:N_mu][i],hypdf[:wq_lo][(i-1)%M+1] )/Float64(T)
    end
	end
	s_typ2err = zeros(Float64,size(hypdf,1))
	for i = 1:size(hypdf,1)
    if hypdf[c_symbol][i] > 0.0  # conformity
      hypdf[:s_typ2err][i] = count_less_s_homo( df, hypdf[c_symbol][i], c_symbol, hypdf[:N][i],hypdf[:N_mu][i],hypdf[:sq_hi][(i-1)%M+1] )/Float64(T)
    elseif hypdf[c_symbol][i] < 0.0  # anti-conformity
      hypdf[:s_typ2err][i] = count_more_s_homo( df, hypdf[c_symbol][i], c_symbol, hypdf[:N][i],hypdf[:N_mu][i],hypdf[:sq_lo][(i-1)%M+1] )/Float64(T)
    end
	end
	sp_t2err = zeros(Float64,size(hypdf,1))
	for i = 1:size(hypdf,1)
    hypdf[:sp_t2err][i] = count_s_prob_accept(df, hypdf[c_symbol][i], c_symbol, hypdf[:N][i], hypdf[:N_mu][i] )/Float64(T)
	end
	#writetable("$(simname)_hypdf.csv",hypdf)
  header_lines = read_headers("$(simname).csv")
  significance_header = " significance level: $lo_val"
  push!(header_lines, significance_header )
  write_dataframe("$(simname)_hypdf.csv",header_lines,hypdf)
	hypdf
end

@doc """ function check_zero_first( N::Int64, hypdf::DataFrame )
Produces an error if the number of "cpower==0.0" rows of the dataframe hypdf is less than N.
""" 
function check_zero_first( N::Int64, hypdf::DataFrame, c_symbol::Symbol )
	i = 1
	while hypdf[c_symbol][i] == 0.0
		i += 1
	end
  if i <= N
		error("hypdf[c_symbol] == 0 rows must come first in hypdf.")
	end
end 


# Call the function hyptest() when the program is invoked from the command line with the "sinname" argument.
if length(ARGS) != 0
  simname = ARGS[1]
  hypdf = hyptest( simname )
	println("hypdf: ",hypdf)
end
