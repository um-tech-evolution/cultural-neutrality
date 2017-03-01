#= The R library "poweRlaw" must be R installed on the machine where this is run.
Author:  David Blasen, May 2016
Revisions: Alden Wright, Oct. to  Dec. 2016
=#
export power_law_estimates, power_law_bootstrap
#using RCall
#using DataFrames
R"library(poweRlaw)"

@doc """ function power_law_estimates(data_vector, threshold_vector)

    Calls the poweRlaw function estimate_xmin to estimate the paramters of data_vector.
    The estimates of the alpha paramter are conditioned the value of xmin.
"""
function power_law_estimates(data_vector, filename::String=""; PNGflag::Bool=true, title::String="", subtitle::String="", top_str::String="" )
    @RCall.rput title
    @RCall.rput subtitle
    xlab_str = "Frequency of variant"
    ylab_str = "Proportion of variants"
    alpha_str = "alpha ="
    xmin_str = "xmin ="
    R"  m_pl =displ$new($data_vector)
        est_pc = estimate_xmin(m_pl)
        est_pl = m_pl$setXmin(est_pc)
        plt = plot(m_pl,main=$title,sub=$subtitle,xlab=$xlab_str,ylab=$ylab_str)
        lines(m_pl,col=2)
        text(5,1.5*min(plt$y),paste($xmin_str,est_pl$xmin))
        text(5,min(plt$y),paste($alpha_str,signif(est_pl$pars,3)))
        mtext($top_str,side=3)
    "
    est_pc = RCall.rcopy(R"est_pc")
    if PNGflag && length(filename) > 0 && filename[end-3:end] == ".png"
      R" 
        png($filename) 
        plot(m_pl,main=$title,sub=$subtitle,xlab=$xlab_str,ylab=$ylab_str)
        lines(m_pl,col=2)
        text(5,5e-3,paste($xmin_str,est_pl$xmin))
        text(5,4e-3,paste($alpha_str,signif(est_pl$pars,3)))
        dev.off()
      "
    end
    #alpha_vals = RCall.rcopy(R"est_scan")
    vals = Dict(
                "xmin" => est_pc[:xmin],
                "alpha" => est_pc[:pars]
                )
    return vals
end

function power_law_bootstrap(data_vector; nthreads::Int64=8, nsims::Int64=5000  )
    R"  m_pl =displ$new($data_vector)
        est_pl = estimate_xmin(m_pl)
        bsp = bootstrap_p( m_pl, xmins=seq(ceiling(est_pl$xmin/2),ceiling(est_pl$xmin),ceiling(est_pl$xmin/4)),no_of_sims=$nsims, threads=$nthreads)
    "
    est_pc = RCall.rcopy(R"est_pl")
    bsp = RCall.rcopy(R"bsp")
    vals = Dict(
                :xmin => est_pc[:xmin],
                :alpha => est_pc[:pars],
                :p => bsp[:p]
                )
    return vals
end

function test_driver()
    # The following reads 2014 US baby names into an data frame.
    bndf = DataFrames.readtable("../data/BabyNames/us_baby_names_by_year2014.csv")
    answer = power_law_estimates(bndf[:GBirths])
    #answer = power_law_estimates(bndf[:GBirths],[262,280,300,320])
    #answer = power_law_estimates([1,2,3,4,5],[1,2,3])
end
