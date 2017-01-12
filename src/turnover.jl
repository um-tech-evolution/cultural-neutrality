using DataStructures

export topKlist, bottomKlist, topKset, bottomKset, turnover, turnover_lists
    add_to_turnover_lists, turnover_lists_to_dataframe

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
  return result[1:K]
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
  return result[1:K]
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
    if Float64(Ylist[i]) > 5.0*N_mu
      deleteat!( Ylist, i )
    end
  end 
  Zlists = Vector{Int64}[]
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
    Float64[]   # list of variances
  ]
end

function add_to_turnover_lists( tl::Vector{Any}, Ylist::Vector{Int64}, Zlists::Vector{Vector{Int64}}, N::Int64,
    mu::Float64 )
  len = length(Ylist)
  N_lst = Int64[ N for i in 1:len ]
  mu_lst = Float64[ mu for i in 1:len ]
  count_lst = [length(Zlists) for i in 1:len ]
  z_lst = zeros(Float64,len)
  sqr_lst = zeros(Float64,len)
  var_lst = zeros(Float64,len)
  for j = 1:length(Zlists)
    for i = 1:len
      z_lst[i] += Zlists[j][i]
      sqr_lst[i] += Zlists[j][i]^2
    end
  end
  for i = 1:len
    mean_lst[i] /= count_lst[i]
    var_lst[i] = (sqr_lst[i]-count_lst[i]*mean_lst[i]^2)/(count_lst[i]-1.0)
  end
  tl[1] = vcat(tl[1],N_lst)
  tl[2] = vcat(tl[2],mu_lst)
  tl[3] = vcat(tl[3],Ylist)
  tl[4] = vcat(tl[4],count_lst)
  tl[5] = vcat(tl[5],z_lst)
  tl[6] = vcat(tl[6],var_lst)
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
    var=tl[6])
end
