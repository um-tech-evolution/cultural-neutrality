# ../experiments/run.jl tests on 3 example configurations:

cd ../experiments
julia -p 2 -L simulation.jl run.jl configs/example1
julia analysis.jl configs/example1

julia -p 2 -L simulation.jl run.jl configs/example2
julia analysis.jl configs/example2

julia -p 2 -L simulation.jl run.jl configs/example3
julia analysis.jl configs/example3
cd ../test
