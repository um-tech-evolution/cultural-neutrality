using DataFrames
using DataStructures
#=
The function bin_logs returns a data frame which can be plotted to generate a power law plot.
The data frame could be written to a CSV file and imported into another program for plotting.
Or the Julia Gadfly package can be used to generate a plot int the form of an SVG image.
For an example of this approach, see the "analysis.jl" program in 
https://github.com/um-tech-evolution/soc-int-experiment.
=#

@doc """ function count_variants( p::Array{Int64,2} )
Returns a data frame with two columns:
:varcounts :  the counts of the number of occurences of each allele that occurs
   in the alleles matrix p.
:logcounts :  the base 10 logartithm of varcounts entries.
Note that the sum the :varcounts column is ngens*N.
"""
function count_variants( plist::PopList )
  counts = counter(Int)
  for pop in plist
    for a in pop
      push!(counts,a)
    end
  end
  vcounts = sort([counts[k] for k in keys(counts)])
  df = DataFrame()
  df[:varcounts] = vcounts
  df[:logcounts] = map(log10,vcounts)
  df
  #=
  ncopies = counter(Int)
  for v in variants
    push!(ncopies,v)
  end
  count_v = [ncopies[k] for k in keys(ncopies)]
  (counts,ncopies,keys(ncopies),count_v)
  (keys(counts),variants)
  #(sort(variants),sort(count_v))
  =#
end

@doc """ function bin_logs( count_vars::DataFrame, nbins::Int64; fix_empty_bins::Bool=false )

count_vars is the dataframe returned by the count_variants function defined above.
Returns bin counts and log10 bin counts of the logcounts.
The resulting data frame can be used to produce a log-log plot of allele counts.
If the neutral model is true, then the allele count data will be fit by a power law and
the log-log plot will be approximately linear.

Returns a data frame with 3 columns:
:lower_bounds :  the lower bounds for each bin
:bins :  the counts the number of logcounts in each bin.
:log_bins :  log10 of the :bins values
Example call:  bin_logs( vardf, 10)
Note:  It may happen that a bin is empty.  Under the default option, the :log_bins value will
be very negative.  If this happens, the suggestion is to run inf_alleles() again, perhaps with
a larger N or ngens.
"""

function bin_logs( count_vars::DataFrame, nbins::Int64; fix_empty_bins::Bool=false )
  bins = zeros(Float64,nbins)
  lower_bounds = zeros(Float64,nbins)
  #max = ceil(10.0*maximum(count_vars[:logcounts]))/10.0
  max = maximum(count_vars[:logcounts])
  println("max: ",max)
  #mult = (nbins-1)/maximum(count_vars[:logcounts])
  mult0 = (nbins-1)/max
  mult1 = max/(nbins-1)
  println("mult0: ",mult0,"  mult1: ",mult1)
  for j = 1:nbins
    lower_bounds[j] = (j-1)*mult1
  end
  for i in 1:size(count_vars)[1]
    l = count_vars[:logcounts][i]
    bin = Int(floor(mult0*(l))+1)
    #println("l: ",l,"  bin: ",bin)
    bins[bin] += 1 
  end
  for j = 1:nbins
    if bins[j] <= 0
      bins[j] = fix_empty_bins ? 1 : eps() # Fixing bins by setting the value to 1 is "lying"
    end
  end
  df = DataFrame()
  df[:lower_bounds] = lower_bounds
  df[:log_bins] = map(log10,bins)
  df[:bins] = bins
  df
end
