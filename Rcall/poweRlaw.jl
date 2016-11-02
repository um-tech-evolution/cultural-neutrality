#= The R library "poweRlaw" must be R installed on the machine where this is run.
Author:  David Blasen, May 2016
=#
export power_law_estimates
using RCall
using DataFrames
R"library(poweRlaw)"

@doc """ function power_law_estimates(data_vector, threshold_vector)

    Calls the poweRlaw function estimate_xmin to estimate the paramters of data_vector.
    The estimates of the alpha paramter are conditioned on the values of xmin in threshold_vector.
"""
function power_law_estimates(data_vector, filename::String="" )

    if length(filename) > 0 && filename[end-3:end] == ".png"
      R" png($filename) "
    end
    R"  m_pl =displ$new($data_vector)
        est_pc = estimate_xmin(m_pl)
        m_pl$setXmin(est_pc)
    "
    est_pc = rcopy(R"est_pc")
    R"
        plot(m_pl)
        lines(m_pl,col=2)
        text(4.0,0.01,paste0("alpha:",est_pc$pars))
        dev.off()
    "
    #alpha_vals = rcopy(R"est_scan")
    vals = Dict(
                "xmin" => est_pc[:xmin],
                "alpha" => est_pc[:pars]
                )
    return vals
end

function test_driver()
    # The following reads 2014 US baby names into an data frame.
    bndf = readtable("../data/BabyNames/us_baby_names_by_year2014.csv")
    answer = power_law_estimates(bndf[:GBirths])
    #answer = power_law_estimates(bndf[:GBirths],[262,280,300,320])
    #answer = power_law_estimates([1,2,3,4,5],[1,2,3])
end
