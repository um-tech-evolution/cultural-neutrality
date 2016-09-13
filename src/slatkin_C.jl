# Julia bindings for slatkin-exact-tools, i. e. Slatkin's C code as refactored by Mark Madsen.
# https://github.com/mmadsen/slatkin-exact-tools
# Requires loading the shared library "slatkin.so".
# For now, be sure to run `make shared` before using.

export SlatkinResult, ewens_montecarlo, SlatkinEnumResult, slatkin_enum

immutable SlatkinResult
  probability::Float64
  theta_estimate::Float64
end

@doc """ function ewens_montecarlo(num_reps::Int32, obs_list::Vector{Int32})
The sampling version of Slatkin "exact" test.
num_reps   is the number of samples.
obs_list    is the given allele count configuration: should be in decreasing order.
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

@doc """ function ewens_montecarlo(num_reps::Int64, obs_list::Config
Convert the arguments from 64 bit integers to 32 bit integers.
"""

function ewens_montecarlo(num_reps::Int64, obs_list::Config)
  return ewens_montecarlo( Int32(num_reps), Int32[ x for x in obs_list ] )
end

immutable SlatkinEnumResult
  prob_slatkin_exact::Float64
  prob_watterson::Float64
  theta_estimate::Float64
end

@doc """ function slatkin_enum( obs_list::Vector{Int64} )
The exact version of Slatkin "exact" test.
obs_list    is the given allele count configuration: should be in decreasing order.
"""
function slatkin_enum( obs_list::Vector{Int64} )
  obs = Int32[ 0; sort(obs_list, rev=true); 0 ]
  
  ccall((:slatkin_enum, "slatkin_enum.so"), SlatkinEnumResult,
      (Ptr{Int32},), obs)
end

@doc """ function slatkin_enum( obs_list::Vector{Int64} )
The exact version of Slatkin "exact" test.
obs_list    is the given allele count configuration: should be in decreasing order.
"""
function slatkin_enum( obs_list::Config )
  obs = Int32[ 0; sort(obs_list, rev=true ); 0 ]
  
  ccall((:slatkin_enum, "slatkin_enum.so"), SlatkinEnumResult,
      (Ptr{Int32},), obs)
end
