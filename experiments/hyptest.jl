#=
Analyzes the CSV file produced by running simulation.jl and run.jl.
This CSV file must include the case where cpower == 0 which is the unbiased null hypothesis.
Tests the null hypothesis that  theta == 0.0 against the alternatives:
   theta < 0.0    for conformist transmission
   theta > 0.0    for anti-conformist transmission
In the resulting table,  wq05 is the 0.05 quantile and wq95 is the 0.95 quantile.
For the null hypothesis  cpower==0, wq05 is critcal theta value for the theta<0 alternative hypothesis,
  and wp95 is the critical theta value for the theta>0 alternative hypothesis.
A type I error is the incorrect rejection of a true null hypothesis, and a type II error is the failure
  to reject a false null hypothesis.
We assume a 5% type I error probability, and compute the type II error probability using the wq05 
  (for conformist) quantile or the wq95 (for anti-conformist) quantile of the null hypothesis distribution.

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

@doc """ function(count_more_w_theta( df, cpower::Float64, N::Int64, N_mu::Float64, wq05::Float64 )

Counts the number of alternative hypothesis outcomes (assuming N, N_mu, cpower as specified)
  which have a w_theta value greater than wq05 where wq05 is the 5% critical value for
  the null hypothesis.
"""
function count_more_w_theta( df, cpower::Float64, N::Int64, N_mu::Float64, wq05::Float64 )
  countnz( (df[:N].==N)&(df[:N_mu].==N_mu)&(df[:cpower].==cpower)&(df[:w_theta].>wq05))
end

@doc """ function count_less_w_theta( df, cpower::Float64, N::Int64, N_mu::Float64, wq95::Float64 )

Counts the number of alternative hypothesis outcomes (assuming N, N_mu, cpower as specified)
  which have a w_theta value less than wq95 where wq05 is the 5% critical value for
  the null hypothesis.
"""
function count_less_w_theta( df, cpower::Float64, N::Int64, N_mu::Float64, wq95::Float64 )
  countnz( (df[:N].==N)&(df[:N_mu].==N_mu)&(df[:cpower].==cpower)&(df[:w_theta].<wq95))
end

@doc """ function(count_more_s_theta( df, cpower::Float64, N::Int64, N_mu::Float64, sq05::Float64 )

Counts the number of alternative hypothesis outcomes (assuming N, N_mu, cpower as specified)
  which have a s_theta value greater than sq05 where sq05 is the 5% critical value for
  the null hypothesis.
"""
function count_more_s_theta( df, cpower::Float64, N::Int64, N_mu::Float64, sq05::Float64 )
  countnz( (df[:N].==N)&(df[:N_mu].==N_mu)&(df[:cpower].==cpower)&(df[:s_theta].>sq05))
end

@doc """ function count_less_s_theta( df, cpower::Float64, N::Int64, N_mu::Float64, sq95::Float64 )

Counts the number of alternative hypothesis outcomes (assuming N, N_mu, cpower as specified)
  which have a s_theta value less than sq95 where sq05 is the 5% critical value for
  the null hypothesis.
"""
function count_less_s_theta( df, cpower::Float64, N::Int64, N_mu::Float64, sq95::Float64 )
  countnz( (df[:N].==N)&(df[:N_mu].==N_mu)&(df[:cpower].==cpower)&(df[:s_theta].<sq95))
end

@doc """ function hyptest_results( df )

Produces a data frame that includes the 5% and 95% quantiles of the null hypothesis distribution
  and the alternative hypothesis distributions.
"""
function hyptest_results( df::DataFrame )
  result_df = by(df, [ :cpower, :N_mu, :N ] ) do d
    DataFrame(
      mean_w_theta=mean(d[:w_theta]),
      wq05 = quantile(d[:w_theta],0.05),
      wq95 = quantile(d[:w_theta],0.95),
      # The following 4 lines can be uncommented to apply to the Slatkin test instead of the Watterson test.
      # Some additional changes are also needed to do this.
      mean_s_theta=mean(d[:s_theta]),
      sq05 = quantile(d[:s_theta],0.05),
      sq05 = quantile(d[:s_theta],0.95),
      #true_theta= d[:true_theta][1],
    )
  end
  #result_df[ df[:cpower].==0.0 ]
  result_df
end

@doc """ function hyptest( simname )

Adds a column "w_typ2err" to the dataframe computed by function hyptest_results().
For the alternative hypothesis rows, this is the type II error corresponding to a type I error
  probability of 5%.
"""
function hyptest( simname )
	include("$(simname).jl")
	M = length(N_list)*length(mu_list)  # number of cases of N times number of cases of N_mu
	df = readtable("$(simname).csv", makefactors=true, allowcomments=true)
	hypdf = hyptest_results( df )
  sort!(hypdf,cols=[order(:cpower,by=abs)])
	check_zero_first( M, hypdf )
	w_typ2err = zeros(Int64,size(hypdf,1))
	for i = 1:size(hypdf,1)
    if hypdf[:cpower][i] > 0.0
      w_typ2err[i] = count_more_w_theta( df, hypdf[:cpower][i], hypdf[:N][i],hypdf[:N_mu][i],hypdf[:wq05][(i-1)%M+1] )
    elseif hypdf[:cpower][i] < 0.0
      w_typ2err[i] = count_less_w_theta( df, hypdf[:cpower][i], hypdf[:N][i],hypdf[:N_mu][i],hypdf[:wq95][(i-1)%M+1] )
    end
	end
	hypdf[:w_typ2err] = w_typ2err
	s_typ2err = zeros(Int64,size(hypdf,1))
	for i = 1:size(hypdf,1)
    if hypdf[:cpower][i] > 0.0
      s_typ2err[i] = count_more_s_theta( df, hypdf[:cpower][i], hypdf[:N][i],hypdf[:N_mu][i],hypdf[:sq05][(i-1)%M+1] )
    elseif hypdf[:cpower][i] < 0.0
      s_typ2err[i] = count_less_s_theta( df, hypdf[:cpower][i], hypdf[:N][i],hypdf[:N_mu][i],hypdf[:sq95][(i-1)%M+1] )
    end
	end
	hypdf[:s_typ2err] = s_typ2err
	writetable("$(simname)_hypdf.csv",hypdf)
	hypdf
end

@doc """ function check_zero_first( N::Int64, hypdf::DataFrame )

Produces an error if the number of "cpower==0.0" rows of the dataframe hypdf is less than N.
""" 
function check_zero_first( N::Int64, hypdf::DataFrame )
	i = 1
	while hypdf[:cpower][i] == 0.0
		i += 1
	end
  if i <= N
		error("cpower == 0 rows must come first in hypdf.")
	end
end 


# Call the function hyptest() when the program is invoked from the command line with the "sinname" argument.
if length(ARGS) != 0
  simname = ARGS[1]
  hypdf = hyptest( simname )
	println("hypdf: ",hypdf)
end
