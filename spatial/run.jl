using SpatialEvolution

function run_trials( simname::AbstractString ) 
  stream = open("$(simname).csv","w")
  println("stream: ",stream)
  sr = SpatialEvolution.spatial_result(N,num_subpops_list[1],num_env_subpops_list[1],ne_list[1],num_attributes, mu, ngens, burn_in,
      horiz_select, circular_variation, extreme_variation_list[1], normal_stddev )
  sr_list_run = SpatialEvolution.spatial_result_type[]
  trial=1
  for num_subpops in num_subpops_list
    for ne in ne_list
      for extreme_variation in extreme_variation_list
        for num_env_subpops in num_env_subpops_list
          if num_env_subpops > 0 && num_env_subpops < maximum(num_subpops_list)
            error("num_env_subpops must be less than the maximum of num_subpops_list")
          end
          sr = SpatialEvolution.spatial_result(N,num_subpops,num_env_subpops,ne,num_attributes, mu, ngens, burn_in,
             horiz_select, circular_variation, extreme_variation, normal_stddev )
          Base.push!(sr_list_run, sr )
          #println("= = = = = = = =")
          #writerow(STDOUT,trial,sr)
        end
      end
    end
  end
  println("===================================")
  sr_list_result = map(spatial_simulation, sr_list_run )
  trial = 1
  writeheader( stream, num_subpops_list, sr )
  writeheader( STDOUT, num_subpops_list, sr )
  for sr_result in sr_list_result
    writerow(stream,trial,sr_result)
    writerow(STDOUT,trial,sr_result)
    trial += 1
  end
end    

if length(ARGS) == 0
  simname = "configs/example2"
else
  simname = ARGS[1]
end
include("$(simname).jl")
println("simname: ",simname)
#println("simtype: ",simtype)
run_trials( simname )
