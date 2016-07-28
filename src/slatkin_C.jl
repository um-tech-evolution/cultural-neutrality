# Julia bindings for slatkin-exact-tools, i. e. Slatkin's C code as refactored by Mark Madsen.
# https://github.com/mmadsen/slatkin-exact-tools
# Requires loading the shared library "slatkin.so".
# For now, be sure to run `make shared` before using.

export SlatkinResult, ewens_montecarlo

immutable SlatkinResult
  probability::Float64
  theta_estimate::Float64
end

@doc """ function ewens_montecarlo(num_reps::Int32, obs_list::Vector{Int32})
The sampling version of Slatkin "exact" test.
num_reps   is the number of samples.
obs_list    is the given allele count configuration.
"""
function ewens_montecarlo(num_reps::Int32, obs_list::Vector{Int32})
  obs = Int32[0; obs_list; 0]

  ccall((:slatkin_mc, "slatkin.so"), SlatkinResult,
      (Int, Ptr{Int32}), num_reps, obs)
end

@doc """ function ewens_montecarlo(num_reps::Int64, obs_list::Vector{Int64})
Convert the arguments from 64 bit integers to 32 bit integers.
"""

function ewens_montecarlo(num_reps::Int64, obs_list::Vector{Int64})
  return ewens_montecarlo( Int32(num_reps), Int32[ x for x in obs_list ] )
end
