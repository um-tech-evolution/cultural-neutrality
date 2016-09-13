#include("../ewens/stewart.jl")
include("slatpart.jl")


CC = gen_configs()  # All allele count configs up to length 8 constructed by hand
PC = cfgs(8,8)      # All allele count configs up to length 8 constructed by computation
C8 = flatten(CC[8]) 
Btbl = BT(8,8)
Btblr = BTr(8,8)
context("CONFIGS") do
  @fact CC --> PC
end

context("Config Probabilities") do
  for i = 1:length(C8)
    @fact Float64(sprob(C8[i])) --> roughly(number_orders(C8[i])*Pcfg(sum(C8[i]),length(C8[i]),C8[i],Btbl)) 
    @fact sprobr(C8[i]) --> number_orders(C8[i])*Pcfg(sum(C8[i]),length(C8[i]),C8[i],Btblr)
  end
end

# May occasionally fail since ewens_montecarlo and slatkin_sample are stochastic
context("Slatkin exact and sample") do
  for i = 1:length(C8)
    sr = ewens_montecarlo( 10000, C8[i] ).probability  # computed by Slatkin's C code
    se = slatkin_exact( C8[i], Btbl )
    ss = slatkin_sample( C8[i], 10000,  Btbl )
    #println("sr: ",sr,"  se: ",se,"  ss: ",ss)
    @fact abs(ss-se) --> less_than( 0.02 )
    @fact abs(sr-ss) --> less_than( 0.02 )
  end
end

# May occasionally fail since watterson_sample is stochastic
context("Watterson exact and sample") do
  for i = 1:length(C8)
    we = watterson_exact( C8[i], Btbl )
    ws = watterson_sample( C8[i], 100000,  Btbl )
    println("we: ",we,"  ws: ",ws)
    @fact abs(we-ws) --> less_than( 0.02 )
  end
end
