# A "sanity test" for innovation.jl and innovation_collection.jl
using Distributions
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

# Tests some of the functionality of innovation_collection.jl under infinite sites
function testall( N::Int64, ngens::Int64, mu::Float64; prob_ext::Float64=0.3, prob_fix::Float64=0.2 )
  ic = innovation_collection(N,1.0)
  i = 1
  g = 1
  g_limit = 1000
  done = false
  while !done && g < g_limit
    println("generation: ",g)
    for index in ic.active
      print("id: ",ic.list[index].identifier," sg:",ic.list[index].start_gen," freq:",ic.list[index].previous_allele_freq)
      println(" sum_counts:",ic.list[index].sum_counts," sum_het:",ic.list[index].sum_heteroz)
    end
    update_innovations!( ic, g, N )
    num_mutations = rand(Binomial(N,mu))
    #println("g: ",g,"  num_mutations: ",num_mutations)
    if g < ngens  # no burn_in for this test
      for j in 1:num_mutations
        fit = 0.50+rand()  # generate fitnesses with a wide range of fitnesses
        println("generating innovation ",i," with fitness: ",fit)
        ic_push!(ic,innovation(i,N,g,fit))
        i += 1
      end
    end
    g += 1
    done = (g > ngens) && length(ic.active) == 0
    if done
      println("final summary")
    end
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
#test(N)
#println("=====================")
srand(3)
testall()
