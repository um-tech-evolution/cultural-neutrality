# A "sanity test" for innovation.jl and innovation_collection.jl
function test(N::Int64)
  id = 1
  ic = innovation_collection(N,1.0)
  # Create 2 innovations on generation 1
  gen = 1
  ic_push!(ic,innovation(id,N,gen))
  id += 1
  ic_push!(ic,innovation(id,N,gen))
  id += 1
  # Updates for generation 2
  gen += 1
  update_innovations!( ic, gen, N )
  # Create innovation on generation 2
  ic_push!(ic,innovation(id,N,gen))
  id += 1
  # Make 1 extinct and 2 fixed on generation 3
  gen = 3
  update_innovations!( ic, gen, N )
  print_summary( ic, print_lists=true )
  return ic
end

# Tests some of the functionality of  innovation.jl and innovation_collection.jl
function testall( N::Int64, ngens::Int64, mu::Float64; prob_ext::Float64=0.3, prob_fix::Float64=0.2 )
  ic = innovation_collection(N,1.0)
  i = 1
  g = 1
  g_limit = 1000
  done = false
  while !done && g < g_limit
    println("generation: ",g)
    update_innovations!( ic, g, N )
    if g <= ngens
      for j = 1:N
        r = rand()
        if r < mu
          println("generating innovation ",i)
          ic_push!(ic,innovation(i,N,g,0.95+0.1*rand()))
          i += 1
        end
      end
    end
    g += 1
    if g == ngens
      println("active: ",ic.active)
    end
    done = (g > ngens) && length(ic.active) == 0
    print_summary( ic, print_lists=true )
  end
  ic
end

function testall()
  testall(N,ngens,mu)
end

include("../src/NeutralCulturalEvolution.jl")

N=5
ngens=4
mu=0.2
test(N)
println("=====================")
testall()
