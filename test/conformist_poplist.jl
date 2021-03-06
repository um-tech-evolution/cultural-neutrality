# Verifed to work on 12/21/16.
using FactCheck
include("../src/NeutralCulturalEvolution.jl")
#using NeutralCulturalEvolution
#include("../src/conformist_poplist.jl")
#include("../src/freq_scaled_fitness.jl")
context("conformist_poplist.jl") do
  N = 10
  ngens = 8
  conformist_prob = 0.2
  anti_conformist_prob = 0.2
  functs = [power_mixed_conformist_poplist, nearly_neutral_power_mixed_conformist_poplist, acerbi_mixed_conformist_poplist]
  context("conformist poplists") do
    for f in functs
      facts("testing $f") do
        N_mu = 2.0
        plist = f(N,N_mu,ngens,conformist_prob,anti_conformist_prob,combine=false)
        @fact length(plist) --> ngens
        plist = f(N,N_mu,ngens,conformist_prob,anti_conformist_prob,combine=true)
        @fact length(plist) --> 1
        @fact length(plist[1]) --> N*ngens
        N_mu = 0.0  # no mutation
        pcnts = pop_counts64(f(N,N_mu,ngens,conformist_prob,anti_conformist_prob,
            combine=true,uniform_start=true,burn_in=0.0)[1])
        @fact length(pcnts) --> 1
      end
    end
  end
end
