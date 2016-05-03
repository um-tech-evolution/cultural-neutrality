using FactCheck

include("../src/infinite_alleles.jl")

context("infinite_alleles.jl") do
  N = 10
  mu = 0.1
  lst_length = 30
  lst = rand( 1:N, lst_length )
  context("function inf_alleles()") do
    ngens = 20
    ntbl = inf_alleles( N, mu, ngens )
    @fact length(ntbl) --> ngens
    @fact length(ntbl[ngens]) --> N
    ctbl = inf_alleles( N, mu, ngens, copy_funct=conformist_copy, K=Int(N/2), C = 0.5 )
    @fact length(ctbl) --> ngens
    @fact length(ctbl[ngens]) --> N
    atbl = inf_alleles( N, mu, ngens, copy_funct=anti_conformist_copy, K=Int(N/2), C = 0.5 )
    @fact length(atbl) --> ngens
    @fact length(atbl[ngens]) --> N
  end
  context("conformist_copy") do
    K = Int(N/2) 
    toplist = topKlist( lst, K )
    @fact length(toplist) -->  K
    for x in toplist
      @fact x in lst --> true
    end
    C = 1.0
    mu = 0.0
    max  = maximum(lst)
    clst = conformist_copy( lst, mu, K, C )
    @fact length(clst) --> lst_length
    for i in 1:length(clst)
      if lst[i] in toplist
        @fact clst[i] --> lst[i]
      end
    end
  end
end
