#!/usr/bin/bash
# The following should run without errors.  

# No checking for correct results is done in this version.
# In fact, the sp_t2err column is not being computed.

cd ../experiments

julia -p 4 -L simulation.jl run.jl configs/example1
julia hyptest.jl configs/example1


julia -p 4 -L simulation.jl run.jl configs/example2
julia hyptest.jl configs/example3


julia -p 4 -L simulation.jl run.jl configs/example3
julia hyptest.jl configs/example3


julia -p 4 -L simulation.jl run.jl configs/example4
julia hyptest.jl configs/example4


julia -p 4 -L simulation.jl run.jl configs/example5
julia hyptest.jl configs/example5


julia -p 4 -L simulation.jl run.jl configs/T100_nn_adv_n50_pm2
julia significance.jl configs/T100_nn_adv_n50_pm2
julia hyptest.jl configs/T100_nn_adv_n50_pm2

cd ../tests
