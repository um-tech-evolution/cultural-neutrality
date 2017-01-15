using DataStructures
using DataFrames

export topKlist, bottomKlist, topKset, bottomKset, turnover, turnover_lists,
    add_to_turnover_lists, turnover_lists_to_dataframe, neutral_turnover_dataframe,
    nearly_neutral_turnover_to_csv, neutral_turnover_to_csv

@doc """ function topKlist( pop::Population, K::Int64 )
Return the list of the K most frequent elements of pop.
If there are more than K of these, break ties arbitrarily
""" 
function topKlist( pop::Population, K::Int64 )
  c = DataStructures.counter(Int64)
  for x in pop
    Base.push!(c,x)
  end
  result = unique(sort( pop, by=x->c[x], rev=true ))
  return result[1:min(K,length(result))]
  #=
  j = K
  max = c[result[K]]
  while j<= length(result) && c[result[j]] == max 
    j += 1
  end
  result[1:j-1]
  =#
end

@doc """ function bottomKlist( pop::Population, K::Int64 )
Return the list of the K least frequent elements of pop.
Break ties arbitrarily.
""" 
function bottomKlist( pop::Population, K::Int64 )
  c = DataStructures.counter(Int64)
  for x in pop
    Base.push!(c,x)
  end
  result = unique(sort( pop, by=x->c[x] ))
  return result[1:min(K,length(result))]
  #=
  If there are more than K of these, return the list of all of them.
  j = K
  min = c[result[K]]
  while j <= length(result) && c[result[j]] == min 
    j += 1
  end
  result[1:j-1]
  =#
end

@doc """ function Kset( pop::Population, K::Int64 )
Return the setof the K most frequent elements of pop.
""" 
function topKset( pop::Population, K::Int64 )
  Set(topKlist( pop, K ))
end

@doc """ function bottomKset( pop::Population, K::Int64 )
Return the setof the K most frequent elements of pop.
""" 
function bottomKset( pop::Population, K::Int64 )
  Set(topKlist( pop, K ))
end

@doc """ function turnover( pop1::Population, pop2::Population, K::Int64 )
The number of alleles entering the toplist plus the number of alleles leaving the toplist
   in the transition from pop1 to pop1.
This is the defintion of Evans and Giametto rather than the definition of Bentley.
Usually, this value is twice Bentley's value.
The exception is when the toplist has less than K elements and an allele leaves without being replaced.
"""
function turnover( pop1::Population, pop2::Population, K::Int64 )
  toplist1 = topKset( pop1, K )
  toplist2 = topKset( pop2, K )
  length(setdiff( toplist1, toplist2 )) + length(setdiff( toplist2, toplist1 ))
end

@doc """ turnover()
Returns the list of turnover values for a sequence poplst of populations where
  the toplist sizes are specified by Ylist which should be in ascending order.
If check_Ylist is true, then remove y values from Ylist that do not satisfy
  the condition of Evans and Giometto for satisfying their formula z = d*mu^a*y^b*N^c
  where a=0.55, b=0.86, c=0.13, d=1.38 assuming neutrality (page 3 of their paper).
"""
function turnover( poplst::Vector{Population}, Ylist::Vector{Int64}, N_mu::Float64, 
      check_Ylist::Bool=true )
  if check_Ylist
    # Remove elements from Ylist that do not satisfy Evans and Giometto's condition that
    #   y < N_mu/0.15 which is implemented conservatively.
    i = length(Ylist)  
    if i > 0 && Float64(Ylist[i]) > 5.0*N_mu
      deleteat!( Ylist, i )
    end
  end 
  Zlists = Vector{Int64}[]
  if length(Ylist) == 0
    return Zlists
  end
  c_prev = DataStructures.counter(Int64)
  for x in poplst[1]
    Base.push!(c_prev,x)
  end
  result_prev = unique(sort( poplst[1], by=x->c_prev[x], rev=true ))
  for i = 2:length(poplst)
    c_next = DataStructures.counter(Int64)
    for x in poplst[i]
      Base.push!(c_next,x)
    end
    result_next = unique(sort( poplst[i], by=x->c_next[x], rev=true ))
    Zlist = zeros(Int64,length(Ylist))
    j = 1
    for y in Ylist
      prev_set = Set(result_prev[1:(min(y,length(result_prev)))])
      next_set = Set(result_next[1:(min(y,length(result_next)))])
      Zlist[j] = length(setdiff( prev_set, next_set )) + length(setdiff( next_set, prev_set ))
      j += 1
    end
    Base.push!(Zlists,Zlist)
    c_prev = c_next
    result_prev = result_next
  end
  Zlists
end

function turnover_lists()
  Any[
    Int64[],    # list of N values
    Float64[],  # list of mu values
    Int64[],    # list of y values
    Int64[],    # list of counts
    Float64[],  # list of means
    Float64[],  # list of variances
    Int64[]     # list of codes for plot colors or symbols
  ]
end

function add_to_turnover_lists( tl::Vector{Any}, Ylist::Vector{Int64}, Zlists::Vector{Vector{Int64}}, N::Int64,
    mu::Float64, code )
  len = length(Ylist)
  N_lst = Int64[ N for i in 1:len ]
  mu_lst = Float64[ mu for i in 1:len ]
  count_lst = [length(Zlists) for i in 1:len ]
  z_lst = zeros(Float64,len)
  sqr_lst = zeros(Float64,len)
  var_lst = zeros(Float64,len)
  code_lst = [code for i in 1:len ]
  for j = 1:length(Zlists)
    for i = 1:len
      z_lst[i] += Zlists[j][i]
      sqr_lst[i] += Zlists[j][i]^2
    end
  end
  for i = 1:len
    z_lst[i] /= count_lst[i]
    var_lst[i] = (sqr_lst[i]-count_lst[i]*z_lst[i]^2)/(count_lst[i]-1.0)
  end
  tl[1] = vcat(tl[1],N_lst)
  tl[2] = vcat(tl[2],mu_lst)
  tl[3] = vcat(tl[3],Ylist)
  tl[4] = vcat(tl[4],count_lst)
  tl[5] = vcat(tl[5],z_lst)
  tl[6] = vcat(tl[6],var_lst)
  tl[7] = vcat(tl[7],code_lst)
  tl
end

function turnover_lists_to_dataframe( tl::Vector{Any} )
  a = 0.55
  b = 0.86
  c = 0.13
  d = 1.38
  DataFrame(
    N=tl[1],
    mu=tl[2],
    y=tl[3],
    x=[d*tl[2][i]^a*tl[3][i]^b*tl[1][i]^c for i = 1:length(tl[1])],
    count=tl[4],
    z=tl[5],
    var=tl[6],
    code=tl[7])
end

# Combines the above functions into one function
function neutral_turnover_dataframe( N_list::Vector{Int64}, mu_list::Vector{Float64}, ngens::Int64, Ylist::Vector{Int64} )
  sort!(Ylist)
  tl = turnover_lists()
  for N in N_list
    println("N: ",N)
    for mu in mu_list
      N_mu = N*mu
      println("N_mu: ",N_mu)
      spl = simple_poplist(N,N_mu,ngens,combine=false)
      #println("length spl: ",length(spl))
      Ylist_tmp = deepcopy(Ylist)
      zl =turnover( spl, Ylist_tmp, N_mu )
      println("length zl: ",length(zl),"  Ylist: ",Ylist_tmp)
      if length(zl[1]) < length(Ylist_tmp)
        Ylist_tmp = Ylist_tmp[1:length(zl[1])]
      end
      add_to_turnover_lists(tl,Ylist_tmp,zl,N,mu)
      println("length tl[1]: ",length(tl[1]))
    end
  end
  turnover_lists_to_dataframe( tl )
end

# Does both nearly neutral and neutral
function nearly_neutral_turnover_dataframe( N_list::Vector{Int64}, mu_list::Vector{Float64}, ngens::Int64, 
    Ylist::Vector{Int64}, dfe::Function )
  sort!(Ylist)
  tl = turnover_lists()
  code = 0    # used to generate colors or plot symbols in plot
  for N in N_list
    println("N: ",N)
    for mu in mu_list
      N_mu = N*mu
      println("N_mu: ",N_mu)
      #  Do nearly neutral
      spl = nearly_neutral_poplist(N,N_mu,ngens,dfe,combine=false)
      Ylist_tmp = deepcopy(Ylist)
      zl =turnover( spl, Ylist_tmp, N_mu )
      if length(zl[1]) < length(Ylist_tmp)
        Ylist_tmp = Ylist_tmp[1:length(zl[1])]
      end
      code += 1
      add_to_turnover_lists(tl,Ylist_tmp,zl,N,mu,code)
      #  Do neutral
      spl = simple_poplist(N,N_mu,ngens,combine=false)
      #println("length spl: ",length(spl))
      Ylist_tmp = deepcopy(Ylist)
      zl =turnover( spl, Ylist_tmp, N_mu )
      println("length zl: ",length(zl),"  Ylist: ",Ylist_tmp)
      if length(zl[1]) < length(Ylist_tmp)
        Ylist_tmp = Ylist_tmp[1:length(zl[1])]
      end
      code += 1
      add_to_turnover_lists(tl,Ylist_tmp,zl,N,mu,code)
    end
  end
  turnover_lists_to_dataframe( tl )
end

function neutral_turnover_to_csv( filename::String, N_list::Vector{Int64}, mu_list::Vector{Float64}, ngens::Int64, Ylist )
  hl = ["neutral","N_list: $(N_list)","mu_list: $(mu_list)","ylist: $(ylist)","ngens: $(ngens)"]
  ndf = neutral_turnover_dataframe(N_list,mu_list,ngens,ylist)
  open(filename,"w") do f
    write_data_frame(f,hl,ndf)
  end
end

# Does both nearly neutral and neutral
function nearly_neutral_turnover_to_csv( filename::String, N_list::Vector{Int64}, mu_list::Vector{Float64}, 
    dfe::Function, ngens::Int64, ylist::Vector{Int64} )
  hl = ["nearly neutral dfe: $(dfe)","N_list: $(N_list)","mu_list: $(mu_list)","ylist: $(ylist)","ngens: $(ngens)"]
  ndf = nearly_neutral_turnover_dataframe(N_list,mu_list,ngens,ylist,dfe)
  open(filename,"w") do f
    write_data_frame(f,hl,ndf)
  end
end
  

