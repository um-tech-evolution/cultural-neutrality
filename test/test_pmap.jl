# Test pmap (parallel map) where the mapped function includes the C function call to Slatkin's test for neutrality. 
# This all works assuming that the LD_LIBRARY_PATH environment variable is set so that slatkin.so is found.
# Run in the directory  cultural-neutrality/test  where  cultural-neutrality corresponds to the GitHub respository
#    of the same name.

@everywhere include("../src/aliases.jl")
@everywhere include("../src/slatkin.jl")
@everywhere include("../src/poplist.jl")
#@everywhere include("../src/NeutralCulturalEvolution.jl")  # Not sure why this doesn't work

function test_pmap(n::Int64)
  clist = Array{Int32,1}[]
  for i = 1:n
    push!(clist,map(x->Int32(x), sort( rand(1:20, 100), rev=true ))) 
  end
  map(x->ewens_montecarlo(Int32(100000),x), clist)
end

function test_pmap( N, mu, ngens )
  pmap(i->ewens_montecarlo(Int32(100000),pop_counts(neutral_poplist(N,mu,ngens)[ngens])),collect(1:10))
end
