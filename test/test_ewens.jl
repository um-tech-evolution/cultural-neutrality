# Test functions from ewens.jl, stewart.jl, slatpart.jl
include("../src/NeutralCulturalEvolution.jl")
using NeutralCulturalEvolution
using FactCheck

facts("Ewens") do
  #include("slatpart.jl")
  include("stewart.jl")
end

