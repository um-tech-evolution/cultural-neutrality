export variant_type, ideal_type
using Distributions
typealias Population Array{Int64,1}
typealias PopList Array{Population,1}

type variant_type
  parent::Int64   # The variant that gave rise to this variant
  innovation::Int64   # The innovation that gave rise to this variant
  fitness::Float64    # The fitness of this variant
  subpop_index::Int64  # index of the containing subpopulation
  attributes::Vector{Float64}   # attributes of the variant
end

type ideal_type
  ideal::Vector{Float64}   # Ideal values for attributes in the environment of this subpop
end

type spatial_result_type
  N::Int64   # meta-population size
  num_subpops::Int64   # number of subpopulations
  num_env_subpops::Int64   # number of "subpopulations" to use for env variation.  0 means use num_subpops.
  ne::Int64  # number of emmigrants in horizontal transfer
  num_attributes::Int64  # number of attributes of a variant
  mu::Float64     # innovation rate
  ngens::Int64  # number of generations after burn-in
  burn_in::Float64
  horiz_select::Bool       # Whether to use selection during horzontal transfer
  circular_variation::Bool    # Whether to vary ideal values in a circular fashion
  extreme_variation::Bool    # Whether to vary ideal values by randomly choosing between high and low values
  normal_stddev::Float64  # standard deviation of normal distribution of mutation perturbations
  fitness_mean::Float64
  fitness_variance::Float64
  attribute_variance::Float64
end

