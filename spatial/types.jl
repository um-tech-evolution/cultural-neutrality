export variant_type, fitness_location_type
using Distributions
typealias Population Array{Int64,1}
typealias PopList Array{Population,1}

type variant_type
  parent::Int64   # The variant that gave rise to this variant
  fitness::Float64    # The fitness of this variant
  fitness_location::Int64  # the fitness location of this variant
  attributes::Vector{Float64}   # attributes of the variant
end

type fitness_location_type
  ideal::Vector{Float64}   # Ideal values for attributes in the environment of this subpop
end

type spatial_result_type
  N::Int64   # meta-population size
  num_subpops::Int64   # number of subpopulations
  num_fit_locations::Int64   # number of "subpopulations" to use for env variation.  0 means use num_subpops.
  ne::Int64  # number of emmigrants in horizontal transfer
  num_attributes::Int64  # number of attributes of a variant
  mu::Float64     # innovation rate
  ngens::Int64  # number of generations after burn-in
  burn_in::Float64
  use_fit_locations::Bool  # Whether to use fitness locations other than subpopulations
  horiz_select::Bool       # Whether to use selection during horzontal transfer
  circular_variation::Bool    # Whether to vary ideal values in a circular fashion
  extreme_variation::Bool    # Whether to vary ideal values by randomly choosing between high and low values
  normal_stddev::Float64  # standard deviation of normal distribution of mutation perturbations
  ideal_max::Float64      # maximum ideal value for circular and extreme ideal values
  ideal_min::Float64      # minimum ideal value for circular and extreme ideal values
  ideal_range::Float64    # range of ideal values for circular and extreme ideal 
  fitness_mean::Float64
  fitness_variance::Float64
  attribute_variance::Float64
end

