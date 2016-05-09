# Julia bindings for slatkin-exact-tools
# For now, be sure to run `make shared` before using.

immutable SlatkinResult
  probability::Cdouble
  theta_estimate::Cdouble
end

function montecarlo(num_reps::Int32, obs_list::Vector{Int32},
    numalleles::Int32, parallel::Int32)
  obs = Int32[0; obs_list; 0]

  ccall((:slatkin_mc, "./slatkin.so"), SlatkinResult,
      (Int, Ptr{Int32}), num_reps, obs)
end

function ewens_montecarlo(num_reps::Int32, obs_list::Vector{Int32})
  parallel::Int32 = 0
  padded_numalleles::Int32 = length(obs_list) + 2

  montecarlo(num_reps, obs_list, padded_numalleles, parallel)
end

#counts1 = Int32[8, 4, 3, 2, 2, 1, 1, 1, 1, 1, 1]
#ewens_montecarlo(Int32(100000), counts1) |> println

#counts2 = Int32[91, 56, 27, 9, 7, 3, 1]
#ewens_montecarlo(Int32(100000), counts2) |> println

