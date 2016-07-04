#=
Generates a CSV file using normal distribtuions that can be used to test hyptest.jl.
Can be run as a command line program, or can be "included" into julia and run by calling 
  function gen_data(simname).
Example command-line:
>  julia gendata.jl configs/gdata
where the file  configs/gdata.jl  includes the paramter settings for the run and 
the file configs/gdata.csv is the generated output CSV file.

See the parameter file  experiments/7_16_16/gdata.jl  for comments on the paramter settings.
=#

using DataFrames

@doc """ function gen_data( ntrials, cpower_list, means_list, stddv, N, N_mu )

See function gen_data(simname) for description of the parameters.
"""

function gen_data( ntrials::Int64, cpower_list::Vector, means_list::Vector, stddv::Float64, N::Int64, N_mu::Float64 )
	df = DataFrame()
	M = ntrials*length(cpower_list)
	cpower = []
	for c in cpower_list
		cpower = vcat( cpower, [c for j = 1:ntrials] )   
	end
	df[:cpower] = cpower
	N_list = [ N for j =1:M]
	df[:N] = N_list 
	N_mu_list = [ N_mu for j = 1:M ]
	df[:N_mu] = N_mu_list
	w_theta = []
	for i = 1:length(cpower_list)
		w_theta = vcat( w_theta, [ stddv*randn() + means_list[i] for j = 1:ntrials ] )
	end
	df[:w_theta] = w_theta
	df
end

@doc """ function gen_data(simname::AbstractString)

simname is the file name (without the ".jl" extension) of the file that contains the parameters.
"""

function gendata(simname) 
	include("$(simname).jl")
  if length(cp_list) != length(mean_list)
    error("The lengths of cp_list and mean_list must be equal in function gen_data.")
  end
  # The values of T, cp_list, mean_list, stddv, N_list, mu_list are specified in  $(simname).jl.
	gen_data( T, cp_list, mean_list, stddv, N_list[1], mu_list[1] )
end

# Handles the command-line invocation.  Note that a command line argument must be supplied.
if length(ARGS) > 0
  simname = ARGS[1]
  df = gendata(simname)
	println("df: ",df)
	writetable("$(simname).csv",df)
end

function read_df(filename)
  df = readtable("$filename.csv", makefactors=true, allowcomments=true)
end
