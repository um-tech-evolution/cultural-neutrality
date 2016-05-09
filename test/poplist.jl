using FactCheck

include("../src/poplist.jl")

context("poplist.jl") do
  pop = [5, 5, 5, 3, 3, 2, 1]
  context("function neutral_poplist()") do
    top2 = topKlist( pop, 2 )
    @fact top2 --> [5, 3]
  end
  context("function popcounts()") do
    pcounts = pop_counts( pop )
    @fact pcounts --> [3,2,1,1]
  end
  context("function neutral_poplist()") do 
    N = 10
    mu = 0.1
    ngens = 20
    plist = neutral_poplist( N, mu, ngens )
    @fact length(plist) --> ngens
    pcounts = pop_counts( plist[1] )
    @fact pcounts --> [ 1 for i=1:N ]
  end
end
