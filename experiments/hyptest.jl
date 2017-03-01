#=
Nearly neutral version:  8/30/16
Mixed conformist version  12/25/16

Analyzes the CSV file produced by running simulation.jl and run.jl, followed by signficance.jl
This CSV file must include the case where cpower == 0 or ?????  TODO:  fix
  which is the unbiased tranmission case.
If the neutral infinite alleles model holds, then the expected homozygosity is 1.0/(1.0+theta)
    where theta = 2*N*mu.
Applies the Watterson and Slatkin homozygosity tests to the null hypothesis that  the homozygosity 
    predicted by the Watterson and Slatkin tests is equal to  1.0/(1.0+theta) where theta = 2*N*mu.
    Note that N and mu are known parameters of the simulation.  
{Watterson,Slatkin} homozygosity == 1/(1+theta) against the alternatives:
   homo > 1/(1+theta)     for conformist transmission
   homo < 1/(1+theta)     for anti-conformist transmission
Also applies the Slatkin "exact" test (more precicely, the Monte Carlo extension) to the
same null and alternative hypotheses.  The Slatkin test returns a p value.

The mean, lo_quantile, hi_quantile, and type 2 error are computed for the following:
  * watterson_homozygosity (w_homoz, wz)
  * slatkin exact test probability (s_prob, sp)
  * slatkin exact probability by my program (se_prob, se) # should be very close to s_prob
  * watterson probability by my program (we_prob, we)  
  * p-homozygosity for the last name in the dataframe (p_?_?)

Note that the w_homoz_sig column of the CSV file produced by significance.jl are 
    1 when w_homoz significantly rejects the null hypothesis  (1-sided test), 
    0 when w_homoz does not significantly reject the null hypothesis  (1-sided test), 
The same is true for the p_?_?_sig column.

The mean_w_sig column is the mean of the w_homoz_sig column of CSV file produced by significance.jl,
and the mean_p_sig column is the mean of the p_?_?_sig column of CSV file produced by significance.jl.
The diff_sig column is the count of how many times these columns are different.

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
include("../src/NeutralCulturalEvolution.jl")
include("dataframe_io.jl")

const list_symbols = [ :dfe_list, :cpower_list, :acpower_list, :topsize_list, :bottomsize_list, :cprob_list, :acprob_list, 
    :acer_flag_list, :bottom_list, :n_list, :N_mu_list ]

@doc """ function count_field_accept( df, hypdf, i::Int64, multi_value_list_syms::Vector{Symbol}, lo_val::Float64, hi_val::Float64, ffield::Symbol )
  Given a cvalue, a value of N, and a value of N_mu, return the fraction of times that the test selected by the ffield.
  If lo_val > 0.0 and hi_val == 1.0, one-sided low.
  If lo_val == 0.0 and hi_val < 1.0, one-sided high.
  If lo_val >0.0 and hi_val < 1.0, two-sided high.
  This is similar to a pivot table analysis in a spreadsheet like Excel
"""
function count_field_accept( df, hypdf, i::Int64, multi_value_list_syms::Vector{Symbol}, lo_val::Float64, hi_val::Float64, ffield::Symbol )
  da = (df[ffield].>lo_val)&(df[ffield].<hi_val)
  for s in multi_value_list_syms
    da &= (df[s] .== hypdf[s][i])
  end
  countnz( da )
end

@doc """ function hyptest_results( df )
Returns a data frame that includes the lo_val and 1.0-lo_val quantiles of the null hypothesis distribution
  and the alternative hypothesis distributions.
Includes zero arrays for the Watterson and Slatkin homozygosity type 2 error and the Slatkin probability type 2 error.
  These will be computed later.
"""
function hyptest_results( df::DataFrame, multi_value_list_syms::Vector{Symbol}, one_value_list_syms::Vector{Symbol},
    lo_val::Float64, hi_val::Float64, p_sig_symbol::Symbol=:none )
  lo_quantile = lo_val > 0.0 ? lo_val : (1.0-hi_val)
  hi_quantile = hi_val < 1.0 ? hi_val : 1.0-lo_val
  println("lo_quantile: ",lo_quantile,"  hi_quantile: ",hi_quantile)
  # Note that c_symbol might be :cpower or :acer_C or :nn_select.  
  result_df = by(df, multi_value_list_syms ) do d
    if p_sig_symbol != :none
      p_symbol = Symbol(string(p_sig_symbol)[1:(end-4)])   # remove "_sig" from the end of p_sig_symbol
    end
    rdf = DataFrame()
    #=
      n = d[:n][1],
      homozyg = 1.0/(1.0 + 2.0*d[:N_mu][1]),
      mean_w_homoz=mean(d[:w_homoz]),
      mean_s_homoz=mean(d[:s_homoz]),
      #wzq_lo = quantile(d[:w_homoz],lo_quantile),
      #wzq_hi = quantile(d[:w_homoz],hi_quantile),
      mean_s_prob=mean(d[:s_prob]),
      #sp_lo = quantile(d[:s_prob],lo_quantile),
      #sp_hi = quantile(d[:s_prob],hi_quantile),
      sp_t2err = 0.0,
      #mean_se_prob=mean(d[:se_prob]),
      #sep_lo = quantile(d[:se_prob],lo_quantile),
      #sep_hi = quantile(d[:se_prob],hi_quantile),
      sep_t2err = 0.0,
      #mean_we_prob=mean(d[:we_prob]),
      #wep_lo = quantile(d[:we_prob],lo_quantile),
      #wep_hi = quantile(d[:we_prob],hi_quantile),
      wep_t2err = 0.0
    )
    =#
    for s in one_value_list_syms
      rdf[ Symbol(string(s)[1:end-5]) ] = eval(s)[1]
    end
    rdf[:mean_K] = mean(d[:K])
    rdf[:mean_w_homoz] = mean(d[:w_homoz])
    rdf[:mean_s_homoz] = mean(d[:s_homoz])
    rdf[:mean_s_prob] = mean(d[:s_prob])
    rdf[:sp_t2err] =  0.0
    #rdf[:sep_t2err] =  0.0
    #rdf[:wep_t2err] =  0.0
    if p_sig_symbol != :none
      rdf[:p_sym_mean]=mean(d[p_symbol])
      rdf[:p_sym_lo] = quantile(d[p_symbol],lo_quantile)
      rdf[:p_sym_hi]= quantile(d[p_symbol],hi_quantile)
      rdf[:p_sym_t2err ]= 0.0
      rdf[:mean_w_sig]=mean(d[:w_homoz_sig])
      rdf[:mean_p_sig]=mean(d[p_sig_symbol])
      rdf[:diff_sig]= countnz( abs(d[p_sig_symbol]-d[:w_homoz_sig]) )
    end
    rdf
    #=
      p_sym_mean=mean(d[p_symbol]),
      p_sym_lo = quantile(d[p_symbol],lo_quantile),
      p_sym_hi = quantile(d[p_symbol],hi_quantile),
      p_sym_t2err = 0.0,
      mean_w_sig=mean(d[:w_homoz_sig]),
      mean_p_sig=mean(d[p_sig_symbol]),
      diff_sig = countnz( abs(d[p_sig_symbol]-d[:w_homoz_sig]) ),
    =#
  end
  result_df
end



@doc """ function hyptest( simname )

Adds columns "w_typ2err" (Watterson) and "s_typ2err" (Slatkin)  to the dataframe computed by function 
  hyptest_results().
For the alternative hypothesis rows, this is the type II error corresponding to a type I error
  probability of lo_val1.
"""
function hyptest( simname )
  global t2error_sides
  global list_symbols
  println("simname: ",simname)
	include("$(simname).jl")
  multi_value_list_syms = Symbol[]
  one_value_list_syms = Symbol[]
  for s in list_symbols
    #println("s: ",s)
    if isdefined(s) 
      #println("  s: ",s)
      if length(eval(s)) > 1
        ssym = Symbol(string(s)[1:end-5])
        push!( multi_value_list_syms, ssym )
      else
        push!( one_value_list_syms, s )
      end
    end
  end
  println("multi_value_list_syms: ",multi_value_list_syms)
  println("one_value_list_syms: ",one_value_list_syms)
  if !isdefined( :t2error_sides )
    t2error_sides = 2     # use 2-sided if not specified.
    println("setting t2error_sides = 2")
  end
  println("t2error_sides: ",t2error_sides)
  if simtype == 1  || simtype == 3  # power conformity
    if minimum( cpower_list ) < 0.0 && maximum( cpower_list ) > 0.0
      error("power conformity parameters must either be all nonnegative or all nonpositive")
    end
    positive_conformity = (maximum( cprob_list ) > 0.0) && ( maximum( cpower_list ) > 0 )   # Boolean variable
    negative_conformity = (maximum( acprob_list ) > 0.0) && ( minimum( acpower_list ) < 0 )   # Boolean variable
  elseif  simtype == 2   # Acerbi conformity
    positive_conformity = (maximum( cprob_list ) > 0.0)    # Boolean variable
    negative_conformity = (maximum( acprob_list ) > 0.0)    # Boolean variable
  else
    positive_conformity = false
    negative_conformity = false
  end
  if t2error_sides == 1  
    if positive_conformity
      lo_val = 0.0
      hi_val = 1.0-sig_level
    elseif negative_conformity
      lo_val = sig_level
      hi_val = 1.0
    else
      error( "t2error_sides should equal 2 for neutral ")
    end
  elseif t2error_sides == 2
    lo_val = sig_level/2.0
    hi_val = 1.0 - lo_val
  end
  println("lo_val: ",lo_val,"  hi_val: ",hi_val)
	M = length(n_list)*length(N_mu_list)  # number of cases of n times number of cases of N_mu
  sigtest = isdefined(:sig_colsyms) && length(sig_colsyms) > 0
  println("sigtest: ",sigtest)
  if sigtest
	  df = readtable("$(simname)_sigtest.csv", makefactors=true, allowcomments=true)
  else
	  df = readtable("$(simname).csv", makefactors=true, allowcomments=true)
  end
  # TODO repace this IF statement
  #=
  if findfirst(names(df),:cpower) > 0   # power model of conformity
    c_symbol = :cpower
  elseif findfirst(names(df),:acer_C) > 0  # Acerbi model of conformity
    c_symbol = :acer_C
  elseif findfirst(names(df),:nn_select) > 0 
    c_symbol = :nn_select 
  else
    error( "file $(simname).csv must contain a column named one of :cpower, :acer_C, or :nn_select ")
  end
  =#
  # TODO:  what is sigtest and p_sig_symbol?
  if sigtest
    p_sig_symbol = names(df)[end]  # assume that the last name is the p_homozygosity symbol
    p_symbol = Symbol(string(p_sig_symbol)[1:(end-4)])   # remove "_sig" from the end of p_sig_symbol
    println("p_sig_symbol: ",names(df)[end])
	  hypdf = hyptest_results( df, multi_value_list_syms, one_value_list_syms, lo_val, hi_val, p_sig_symbol )
  else
	  hypdf = hyptest_results( df, multi_value_list_syms, one_value_list_syms, lo_val, hi_val )
  end
  # sort hypdf on all multi_value columns in order given by list_symbols at the beginning of this file.
  sort_sym_list = [ (s==:dfe ? order(s, lt=(x,y)->(string(x) < string(y))): order(s, by=abs)) for s in multi_value_list_syms]
  sort!(hypdf,cols=sort_sym_list)  
  #println("hypdf:")
  #println(hypdf)
	#check_zero_first( M, hypdf, c_symbol )
  if sigtest
    #fields = [:s_prob, :se_prob, :we_prob, p_symbol]
    #t2_fields = [:sp_t2err, :sep_t2err, :wep_t2err, :p_sym_t2err ]  #TODO: define these automatically
    fields = [:s_prob, p_symbol]
    t2_fields = [:sp_t2err, :p_sym_t2err ]  #TODO: define these automatically
  else
    #fields = [:s_prob, :se_prob, :we_prob ]
    #t2_fields = [:sp_t2err, :sep_t2err, :wep_t2err ]  #TODO: define these automatically
    fields = [:s_prob]
    t2_fields = [:sp_t2err]
  end
  for j = 1:length(fields)
	  for i = 1:size(hypdf,1)
      #hypdf[t2_fields[j]][i] = count_field_accept(df, lo_val, hi_val, fields[j], hypdf[c_symbol][i], c_symbol, hypdf[:n][i], hypdf[:N_mu][i] )/Float64(T)
      hypdf[t2_fields[j]][i] = count_field_accept(df, hypdf, i, multi_value_list_syms, lo_val, hi_val, fields[j]  )/Float64(T)
      #hypdf[t2_fields[j]][i] = count_field_accept(df, lo_2sided, hi_2sided, fields[j], hypdf[c_symbol][i], c_symbol, hypdf[:n][i], hypdf[:N_mu][i] )/Float64(T)
    end
	end
  if sigtest
    # Rename the p-homozygsity columns that correspond to p_sig_symbol (see above) to be more meaningful.
    p_sym_names = [:p_sym_mean, :p_sym_lo, :p_sym_hi, :p_sym_t2err]
    p_symbol = Symbol(string(p_sig_symbol)[1:(end-4)])   # remove "_sig" from the end of p_sig_symbol
    p_sym_string = string(p_symbol)
    p_new_names = [ Symbol("$(p_sym_string)_mean"), Symbol("$(p_sym_string)_lo"), Symbol("$(p_sym_string)_hi"),
      Symbol("$(p_sym_string)_t2_err") ]
    rename!( hypdf, p_sym_names, p_new_names )
  end
	#writetable("$(simname)_hypdf.csv",hypdf)
  if sigtest
    header_lines = read_headers("$(simname)_sigtest.csv")
  else
    header_lines = read_headers("$(simname).csv")
  end
  write_dataframe("$(simname)_hyptest.csv",header_lines,hypdf)
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
