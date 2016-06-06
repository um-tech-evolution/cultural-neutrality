#module NeutralCulturalAnalysis

using DataFrames
using Gadfly

if length(ARGS) == 0
  error("No simulation config specified.")
end

simname = ARGS[1]

include("$(simname).jl")


function slatkin_prob_results( df )
  result_df = by(df, [:N_mu, :N ] ) do d
    DataFrame(
      mean_prob=mean(d[:prob]),
      median_prob=median(d[:prob]),
      q01=findlast(x->x<0.01,sort(d[:prob])),
      q025=findlast(x->x<0.025,sort(d[:prob])),
      q05=findlast(x->x<0.05,sort(d[:prob])),
      q95=T-findfirst(x->x>0.95,sort(d[:prob])),
      q975=T-findfirst(x->x>0.975,sort(d[:prob])),
      q99=T-findfirst(x->x>0.99,sort(d[:prob])),
    )
  end
  result_df
end

function theta_results( df )
  result_df = by(df, [:N_mu, :N ] ) do d
    DataFrame(
      mean_theta=mean(d[:est_theta]),
      median_theta=median(d[:est_theta]),
      true_theta= d[:true_theta][1],
      #q01 = quantile(d[:est_theta],0.01),
      #q025= quantile(d[:est_theta],0.025),
      #q05 = quantile(d[:est_theta],0.05),
      #q95 = quantile(d[:est_theta],0.95),
      #q975 = quantile(d[:est_theta],0.975),
      #q99 = quantile(d[:est_theta],0.99)
    )
  end
  result_df
end

function K_results( df )
  result_df = by(df, [:N_mu, :N ] ) do d
    DataFrame(
      mean_K_est=mean(d[:K]),
      median_K_est=median(d[:K]),
      true_theta= d[:true_theta][1],
      #q01 = quantile(d[:est_theta],0.01),
      #q025= quantile(d[:est_theta],0.025),
      #q05 = quantile(d[:est_theta],0.05),
      #q95 = quantile(d[:est_theta],0.95),
      #q975 = quantile(d[:est_theta],0.975),
      #q99 = quantile(d[:est_theta],0.99)
    )
  end
  result_df
end


# Plot development

function drawPlot(variable)
  for kval = collect(K)
    fdf = df[df[:K] .== kval,:]
    println("kval: ",kval,"  fdf: ",fdf)
    p = plot(fdf, x="generation", y=variable, color="siN_mutationType",
      Geom.line,
      Guide.title("$(variable) (N=$(N), K=$(kval), SEL=$(S), Trials=$(T))"))
    kvalstr = @sprintf("%02d", kval)
    draw(SVG("$(simname)_$(variable)_K$(kvalstr).svg", 8inch, 8inch), p)
  end
end


df = readtable("$(simname).csv", makefactors=true, allowcomments=true)
#println("df: ",df)
if simtype == 1
  spdf = slatkin_prob_results( df )
  println("spdf: ",spdf)
  writetable("$(simname)_slat_prob.csv",spdf)
  stdf = theta_results( df )
  println("thetadf slatkin: ",stdf)
  writetable("$(simname)_slatkin_theta.csv",stdf)
elseif simtype == 2  # Watterson
  stdf = theta_results( df )
  println("thetadf watterson: ",stdf)
  writetable("$(simname)_watterson_theta.csv",stdf)
elseif simtype == 3 # K estimate
  kdf = K_results( df )
  println("Kdf  ",kdf)
  writetable("$(simname)_K.csv",kdf)
end


#=
drawPlot("meanFitness")
drawPlot("medianFitness")
drawPlot("maxFitness")
=#

#end

