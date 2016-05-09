include("../src/conformist_poplist.jl")
include("../src/freq_scaled_fitness.jl")
context("conformist_poplist.jl") do
  N = 10
  mu = 0.1
  ngens = 8
  context("power_conformist_poplist()") do
    plist = power_confirmist_poplist(N,mu,ngens,conformist_power=1.0)
    @fact length(plist) --> ngens
  end
  context("function acerbi_conformist_poplist()") do
    # TOTO
  end
end
