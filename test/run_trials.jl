using FactCheck
include("../src/NeutralCulturalEvolution.jl")
#using NeutralCulturalEvolution
#include("../src/conformist_poplist.jl")
#include("../src/freq_scaled_fitness.jl")
include("../experiments/mixed_conformist.jl")
context("run_trials.jl") do
  n = N = 250; N_mu = 2.0; ngens = 1000; burn_in = 2.0;
  context("run_trials_power_mixed_conformist()") do
    cprob=0.2; acprob=0.2;
    cpwr = 0.5; acpwr = -0.5;
    pmv = run_trials_power_mixed_conformist(n,N,N_mu,ngens,cprob,acprob,conformist_power=cpwr,anti_conformist_power=acpwr,burn_in=burn_in,CSVflag=true)
    println(pmv)
# Run nearly neutral with default settings
    #@fact length(plist) --> ngens
  end
  context("function run_trials_acerbi_mixed_conformist_poplist()") do
    tsize=3; bsize=1;
    amc1 = run_trials_acerbi_mixed_conformist(n,N,N_mu,ngens,cprob,acprob,toplist_size=tsize,bottomlist_size=bsize,bottom=false,burn_in=burn_in,CSVflag=true)
    println(amc1)
    amc2 = run_trials_acerbi_mixed_conformist(n,N,N_mu,ngens,cprob,acprob,toplist_size=tsize,bottomlist_size=bsize,bottom=true,burn_in=burn_in,CSVflag=true)
    println("amc2[alpha]: ",amc2["alpha"])
    println("amc2[xmin]: ",amc2["xmin"])
  end
  context("function run_trials_nearly_neutral()") do
    rtnn_adv = run_trials_nearly_neutral(n,N,N_mu,ngens,dfe=dfe_advantageous)
    rtnn_del = run_trials_nearly_neutral(n,N,N_mu,ngens,dfe=dfe_deleterious)
    rtnn_mix = run_trials_nearly_neutral(n,N,N_mu,ngens,dfe=dfe_mixed)
    rtnn_mod = run_trials_nearly_neutral(n,N,N_mu,ngens,dfe=dfe_mod)
  end
  context("function run_trials_nearly_neutral_power_mixed_conformist()") do
    rtnn_pmc = run_trials_nearly_neutral_power_mixed_conformist(n,N,N_mu,ngens,cprob,acprob,conformist_power=cpwr,anti_conformist_power=acpwr,burn_in=burn_in,CSVflag=true,dfe=dfe_advantageous,alpha=1.0,theta=0.01) 
  end
end
