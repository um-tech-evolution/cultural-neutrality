using DataFrames
export cumm_dist, significance_point
#include("../experiments/hyptest.jl")
include("../src/NeutralCulturalEvolution.jl")
include("dataframe_io.jl")
#import ProgressMeter
#const PM = ProgressMeter


@doc """ function cumm_dist( clist::ConfigList, Btbl::Array{Float64,2}, homoz_funct::Function; normalize::Bool=false )
Starting with a list clist of allele count configurations (as produced by function ncfgs()), 
    compute the Stewart probability of each configuraiton (by applying function Pcfg), and compute
    the p-homozygosity (by applying function homoz_funct), and then sort by p-homozygosity.  
    Convert to cummulative probabilities, and return a cummulative probability list of triples, where the elements 
    of the triple are:
   1)  the allele count configuration
   2)  the cummulative probability of the configuration.
   3)  the p-homozygosity of the configuration.
If the flag normalize == true, then normalize so that the last cummulative probability is 1.0.
This is helpful if clist was produced by simulation and low probability configurations are under-represented.
"""
function cumm_dist( clist::ConfigList, Btbl::Array{Float64,2}, homoz_funct::Function; normalize::Bool=false )
  N = convert(Int64,sum(clist[1]))
  K = length(clist[1])
  cummP = map(c->(c,number_ordersF(c)*Pcfg(N,K,c,Btbl),homoz_funct(c)),clist)   # table of (config,probability,homoz) tuples
  sort!(cummP,by=c->c[3])
  for i in 2:length(cummP)
    cummP[i] = (cummP[i][1],cummP[i][2]+cummP[i-1][2],cummP[i][3])
  end
  if normalize &&  cummP[end][2] < 1.0-10.0*eps()
    multiplier = 1.0/cummP[end][2]
    println("normalizing multiplier: ",multiplier)
    for i in 1:length(cummP)
      cummP[i] = (cummP[i][1],multiplier*cummP[i][2],cummP[i][3])
    end
  end
  cummP
end

@doc """ function significance_point( sig_cutoff::Float64, cummP::Array{Tuple{Config,Float64,Float64},1} )
Returns a "signfificance point" for the p-homozygosity statistic used to construct the cummulative probability array  cummP.
sig_cutoff  should be a number between 0.0 and 1.0.  
If sig_cutoff==0.0, the minimum p-homozygosity is returned.
If sig_cutoff==1.0, the maximum p-homozygosity is returned.
If sig_cutoff==s with 0 < s <= 0.5, then the minimum p-homozygosity q  such that Pr( h <= q ) == s is returned.
If sig_cutoff==s with 0.5 < s < 1.0, then the maximum p-homozygosity q  such that Pr( q <= h ) == s is returned.
These probabilities on h are computed assuming the neutral infinite alleles model with sample size N and allele count K.
"""
function significance_point( sig_cutoff::Float64, cummP::Array{Tuple{Config,Float64,Float64},1} )
  if sig_cutoff == 0.0
    return cummP[1][3]
  elseif sig_cutoff == 1.0
    return cummP[end][3]
  elseif sig_cutoff <= 0.5
    j = findlast(x->x[2]<sig_cutoff,cummP)
    if j == 0
      return "-"
    else
      return cummP[j][3]
    end
  else
    j = findfirst(x->x[2]>sig_cutoff,cummP)
    if j == 0 || j == length(cummP)
      return "-"
    else
      return cummP[j][3]
    end
  end
end

@doc """ function significance_dicts( sig_cutoff::Float64, n::Int64, funct_list::Vector{Function} )
Returns a list of dictionaries, one for each function in funct_list, of sig_cutoff sigificance points indexed on k.
In other words, dict is the returned list, dict[i][k] is the significance point for function funct_list[i]
corresponding to number of alleles k.
See the comments on function significance_point for the defintion of the significance point corresponding to sig_cutoff.
"""
function significance_dicts( sig_cutoff::Float64, n::Int64, funct_list::Vector{Function} )
  #progress = PM.Progress(T, 1, "Running...", 40)
  Btbl = BT( n, n )
  dict_list = Dict{Int64,Float64}[]
  for i = 1:length(funct_list)
    sig_dict = Dict{Int64,Float64}()
    push!(dict_list,sig_dict)
  end
  for k = 1:n
    C = ncfgs(n,k)
    #println("significance_dicts: k: ",k)
    for i = 1:length(funct_list)
      cP =     cumm_dist( C, Btbl, funct_list[i] )
      sp = significance_point( sig_cutoff, cP )
      if sp == "-"
        dict_list[i][k] = -1.0
      else
        dict_list[i][k] = sp
      end
    end
  end
  dict_list
end

function build_allele_configs( n::Int64 )
  #configs_list = Vector{ConfigList}()
  configs_list = pmap(k->ncfgs(n,k), collect(1:n))
end

function significance_dicts_map( sig_cutoff::Float64, n::Int64, funct_list::Vector{Function} )
  Btbl = BT( n, n )
  dict_list = Dict{Int64,Float64}[]
  for i = 1:length(funct_list)
    sig_dict = Dict{Int64,Float64}()
    push!(dict_list,sig_dict)
  end
  Clist = build_allele_configs(n)
  for k = 1:n
    #println("significance_dicts: k: ",k)
    for i = 1:length(funct_list)
      cP =     cumm_dist( Clist[k], Btbl, funct_list[i] )
      sp = significance_point( sig_cutoff, cP )
      if sp == "-"
        dict_list[i][k] = -1.0
      else
        dict_list[i][k] = sp
      end
    end
  end
  dict_list
end

@doc """ function add_significance_columns( df::DataFrame,  sig_cutoff::Float64, n::Int64, col_symbol_list::Vector{Symbol}, 
    funct_list::Vector{Function}, dict_list::Vector{Dict{Int64,Float64}} )
  Adds 0/1 columns corresponding to Watterson homozygosity and/or p-homozygosity.
     A 1 value shows significance, and 0 value shows no signficance.
     As of 8/23/16, only one-sided signficance is impllemented.
"""
function add_significance_columns( df::DataFrame,  sig_cutoff::Float64, n::Int64, col_symbol_list::Vector{Symbol}, 
    funct_list::Vector{Function}, dict_list::Vector{Dict{Int64,Float64}} )
  if length(col_symbol_list) != length(funct_list)
    error("col_symbol_list and funct_list must be of the same length in function add_significance_columns")
  end
  #dict_list = significance_dicts( sig_cutoff, n, funct_list )
  for i = 1:length(col_symbol_list)
    #println("i: ",i,"  col_sym: ",col_symbol_list[i])
    new_column = zeros(Int64,length(df[col_symbol_list[i]]))
    for j = 1:length(df[col_symbol_list[i]])
      k = df[:K][j]
      if sig_cutoff > 0.5
        new_column[j] =  df[col_symbol_list[i]][j]>dict_list[i][k] ? 1 : 0 
      else
        new_column[j] =  df[col_symbol_list[i]][j]<dict_list[i][k] ? 1 : 0 
      end
    end
    df[Symbol(col_symbol_list[i],"_sig")] = new_column
    #println("i:",i,"  new_column: ",transpose(new_column))
  end
  df
end

@doc """ function write_sig_dataframe( simname )
 Reads the dataframe produced by the simulation, then adds signficance columns to this dataframe,
    and writes the dataframe to a file.
 The parameters are specified in the file  "\$(simname).jl".
 The sig_level parameter should be less than 0.5.
"""
function write_sig_dataframe( simname )
  include("$(simname).jl")
  sig_functs = map(eval,sig_colsyms)   # eval'ing a symbol returns the value, which in this case is a function
  df = readtable("$(simname).csv", makefactors=true, allowcomments=true)
  if minimum(df[:n]) != maximum(df[:n])
    error("Dataframe must have an :n column with a single value.")
  end
  n = minimum(df[:n])
  println("sample size n:",n)
  if isdefined(:cp_list)  # power conformity
    if minimum( cp_list ) < 0.0 && maximum( cp_list ) > 0.0
      error("power conformity parameters must either be all nonnegative or all nonpositive")
    end
    if maximum( cp_list ) > 0  # positive power conformity
      sig_cutoff = 1.0 - sig_level
    elseif  minimum(cp_list) < 0  # negative power conformity
      sig_cutoff = sig_level
    end
  elseif isdefined(:acer_C_list)  # Acerbi conformity
    if minimum( acer_C_list ) < 0.0 && maximum( cp_list ) > 0.0
      error("Acerbi conformity parameters must either be all nonnegative or all nonpositive")
    end
    if maximum( acer_C_list ) > 0  # positive power conformity
      sig_cutoff = 1.0 - sig_level
    elseif  minimum(acer_C_list) < 0  # negative power conformity
      sig_cutoff = sig_level
    end
  elseif isdefined(:dfe)  # Nearly Neutral
    sig_cutoff = sig_level
  end
  sig_dict_list = significance_dicts(sig_cutoff,n,sig_functs)
  #println(sig_dict_list[1])
  sig_df = add_significance_columns( df, sig_cutoff, n, sig_colsyms, sig_functs, sig_dict_list)
  header_lines = read_headers("$(simname).csv")
  significance_level_header = " significance level: $sig_level"
  push!(header_lines, significance_level_header )
  colsyms_header = " sig_colsyms: $(sig_colsyms)"
  push!(header_lines, colsyms_header )
  functs_header = " sig_functs: $(sig_functs)"
  push!(header_lines, functs_header )
  significance_cutoff_header = " significance cutoff: $sig_cutoff"
  push!(header_lines, significance_cutoff_header )
  write_dataframe("$(simname)_sigtest.csv",header_lines,sig_df)
end


# Call the function write_sig_dataframe() when the program is invoked from the command line with the "simname" argument.
if length(ARGS) != 0
  simname = ARGS[1]
  println("simname: ",simname)
  write_sig_dataframe( simname )
  #println("hypdf: ",hypdf)
end
