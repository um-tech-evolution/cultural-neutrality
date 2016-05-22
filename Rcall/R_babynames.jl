# The R library "poweRlaw" must be R installed on the machine where this is run.
using RCall
using DataFrames
R"library(poweRlaw)"

@doc """ function power_law_estimates(data_vector, threshold_vector)

    Calls the poweRlaw function estimate_xmin to estimate the paramters of data_vector.
    The estimates of the alpha paramter are conditioned on the values of xmin in threshold_vector.
"""
function power_law_estimates(data_vector, threshold_vector)
    R"m_pl =displ$new($data_vector)"
    xmin = rcopy(R"estimate_xmin(m_pl)")[:xmin]
    R"
        threshold_array = $threshold_vector
        est_scan = 0*threshold_array
        for(i in seq_along(threshold_array)){
            m_pl$setXmin(threshold_array[i])
            est_scan[i] = estimate_pars(m_pl)$pars
        }
    "
    alpha_vals = rcopy(R"est_scan")
    vals = Dict(
                "xmin" => xmin,
                "alpha_vector" => alpha_vals
                )
    return vals
end

function test_driver()
    # The following reads 2014 US baby names into an data frame.
    bndf = readtable("../data/BabyNames/us_baby_names_by_year2014.csv")
    answer = power_law_estimates(bndf[:GBirths],[262,280,300,320])
    #answer = power_law_estimates([1,2,3,4,5],[1,2,3])
end
