function _copy(G::Graph{Undirected})
  GG = Graph{Undirected}(n_vertices(G))
  for e in edges(G)
    add_edge!(GG, src(e), dst(e))
  end
  return GG
end

# multi-graph type for Kontsevich moduli spaces
struct MultiGraph
  g::Graph{Undirected}             # underlying SimpleGraph
  multi::Dict{Edge, Int} # multiplicity of each edge as a tuple (a,b)
  aut::Int                      # number of automorphisms
end

Base.show(io::IO, g::MultiGraph) = print(io, "MultiGraph based on ", g.g)

function Base.getindex(g::MultiGraph, e::Edge)
  if haskey(g.multi, e)
    return g.multi[e]
  else haskey(g.multi, Edge(dst(e), src(e)))
    return g.multi[Edge(dst(e), src(e))]
  end
  g.multi[e]
  #s, d = src(e), dst(e)
  #g.multi[(min(s, d), max(s, d))]
end

Oscar.edges(g::MultiGraph) = edges(g.g)
Oscar.vertices(g::MultiGraph) = vertices(g.g)
Oscar.all_neighbors(g::MultiGraph, v) = all_neighbors(g.g, v)

# type of a variety (or orbifold) with a torus action
###@attributes mutable struct TnVariety{P}
 ### dim::Int
 ### points::Vector{Pair{P, Int}} # pairs of (point => orbifold multiplicity)

 ### function TnVariety(n::Int, points::Vector{Pair{P, Int}}) where P
 ###   new{P}(n, points)
 ### end
###end

###Base.show(io::IO, X::TnVariety) = print(io,
###  "TnVariety of dim ", X.dim, " with ", length(X.points), " fixed points")

struct Cycle
  loc::Function
  num_alloc::Int
end

(c::Cycle)(args...) = c.loc(args...)

function Base.:*(c::Cycle, d::Cycle)
  n = max(c.num_alloc, d.num_alloc) + 1
  loc(p::Tuple{MultiGraph, Vector}, λ::Vector,
    ans::QQFieldElem=QQ(1), alloc::Vector{QQ}=[QQ() for _ in n]) = begin
    Nemo.set!(ans, 1, 1)
    c.loc(p, λ, ans, alloc)
    d.loc(p, λ, alloc[end], alloc)
    return muleq!(ans, alloc[end])
  end
  Cycle(loc, n)
end

###############################################################################
# 
# main interface for Kontsevich moduli space
# 

@doc raw"""
     kontsevich_moduli_space(n::Int, d::Int; weights=nothing)

Return the Kontsevich moduli space $\overline{M}(\mathbb P^n, d)$ as a `TnVariety` (abstract orbifold).

!!! note
    With respect to notation, we refer to [Dan14](@cite).

# Examples
```jldoctest
julia> K = kontsevich_moduli_space(4, 3)
TnVariety of dim 16 with 740 fixed points

julia> V = fixed_points(K);

julia> mult = V[740][2]
6

```
"""
function kontsevich_moduli_space(n::Int, d::Int; weights=nothing)
  N = (n+1)*(d+1) - 4
  points = [(g, c) => g.aut for g in multi_trees(d) for c in colors(g.g, n)]
  M = TnVariety(N, points)
  set_attribute!(M, :dim => n)
  #set_special(M, :dim => n)
  set_attribute!(M, :deg => d)
  #set_special(M, :deg => d)
  set_attribute!(M, :weights => weights)
  #set_special(M, :weights => weights)
  M
end

function integral(X::TnVariety, cycle::Cycle; noretry::Bool=false)
  n = get_attribute(X, :dim)
  #n = get_special(X, :dim)
  λ = get_attribute(X, :weights)
  #λ = get_special(X, :weights)
  if λ === nothing
    λ = Int[rand(Int16) for _ in 0:n]
  end
  retry = true
  while retry
    if noretry
      retry = false
    end
    # pre-allocations
    # this is wrong and needs to be fixed if we want a parallel version
    # _ans = [QQ() for _ in 1:Threads.nthreads()]
    # extra = [(QQ(), QQ()) for _ in 1:Threads.nthreads()]
    # alloc = [[QQ() for _ in 1:max(cycle.num_alloc, 3)] for _ in 1:Threads.nthreads()]
    _ans = QQ()
    extra = (QQ(), QQ())
    alloc = [QQ() for _ in 1:max(cycle.num_alloc, 3)]
    try
      # the main loop using Bott's formula
      #Threads.@threads
      for (p, e) in X.points
        Fp, Tp = extra
        cycle(p, λ, Fp, alloc)
        __euler(p, λ, Tp, alloc)
        Nemo.add!(_ans, div!(Fp, mul!(Tp, e)))
      end
    catch e
      if e isa InterruptException
        rethrow(e)
      end
      !noretry || throw(e)
      # It's possible that the chosen random weights do not work
      # In this case we reset and start again
      for i in 1:n+1
        λ[i] = rand(Int16)
      end
      continue
    end
    ans = QQ()
    #for x in _ans
    #  Oscar.addeq!(ans, x)
    #end
    return _ans
  end
end

hypersurface(ns::Int...) = hypersurface(collect(ns))

function hypersurface(ns::Vector{Int})
  loc(p::Tuple{MultiGraph, Vector}, λ::Vector,
    ans::QQFieldElem=QQ(1), alloc::Vector{QQFieldElem}=[QQ(), QQ()]) = begin
    Nemo.set!(ans, 1, 1)
    a1, a2 = alloc
    t = _tally(ns)
    for n in keys(t)
      _hypersurface(p, λ, n, a1, [a2])
      mul!(ans, pow!(a1, t[n]))
    end
    return ans
  end
  Cycle(loc, 2)
end

function _hypersurface(p::Tuple{MultiGraph, Vector}, λ::Vector, n::Int,
    ans::QQFieldElem=QQ(1), alloc::Vector{QQFieldElem}=[QQ()])
  g, i = p
  Nemo.set!(ans, 1, 1)
  a1, = alloc
  for e in edges(g)
    for a in 0:n*g[e]
      b = n*g[e] - a
      Nemo.mul!(ans, Nemo.set!(a1, a*λ[i[src(e)]] + b*λ[i[dst(e)]], g[e]))
    end
  end
  for v in vertices(g)
    Nemo.set!(a1, n*λ[i[v]], 1)
    pow!(a1, 1 - length(all_neighbors(g, v)))
    mul!(ans, a1)
  end
  ans
end

function _euler(p::Tuple{MultiGraph, Vector}, λ::Vector, ans::QQFieldElem=QQ(1),
    alloc::Vector{QQFieldElem}=[QQ(), QQ(), QQ()])
  g, i = p
  Nemo.set!(ans, 1, 1)
  n = length(λ) - 1
  a1, a2, a3 = alloc

  for e in edges(g)
    Nemo.set!(a1, λ[i[src(e)]] - λ[i[dst(e)]], g[e])
    pow!(a1, 2*g[e])
    mul!(a1, (-1)^g[e] * factorial(g[e])^2)
    mul!(ans, a1)
    for a in 0:g[e]
      b = g[e]-a
      for k in 1:n+1
        if k != i[src(e)] && k != i[dst(e)]
          Nemo.set!(a1, a*λ[i[src(e)]] + b*λ[i[dst(e)]], g[e])
          mul!(ans, sub!(a1, λ[k]))
        end
      end
    end
  end
  for v in vertices(g)
    nbrs = all_neighbors(g, v)
    Nemo.set!(a1, 1, 1)
    for j in 1:n+1
      if j != i[v]
        mul!(a1, λ[i[v]] - λ[j])
      end
    end
    mul!(ans, pow!(a1, 1 - length(nbrs)))
    Nemo.set!(a1, 0, 1)
    Nemo.set!(a2, 1, 1)
    for w in nbrs
      e = Edge(v, w)
      Nemo.set!(a3, λ[i[v]] - λ[i[w]], g[e])
      mul!(a2, a3)
      add!(a1, Nemo.inv!(a3))
    end
    mul!(a2, pow!(a1, 3 - length(nbrs)))
    mul!(ans, a2)
  end
  ans
end

const __euler = Cycle(_euler, 3)

###############################################################################
# 
# functions handling graphs
# 

# number of unlabeled rooted trees https://oeis.org/A000081
function oeis_A81(n::Int)
  n == 1 && return 1
  ans = 0
  for p in artitions(n-1)
    ans += prod(binomial(oeis_A81(k)+c-1, c) for (k,c) in _tally(p))
  end
  ans
end

# number of unlabeled trees https://oeis.org/A000055
function oeis_A55(n::Int)
  n == 1 && return 1
  ans = iseven(n) ? binomial(oeis_A81(n÷2)+1, 2) : 0
  for p in _partitions(n-1, n-1, (n-1)÷2)
    ans += prod(binomial(oeis_A81(k)+c-1, c) for (k,c) in _tally(p))
  end
  ans
end

# unlabeled rooted trees
function oeis_A81_graphs(n::Int)
  n == 1 && return [Graph{Undirected}(1)]
  ans = Graph{Undirected}[]
  for p in partitions(n-1)
    for g in Base.Iterators.product([_sym(oeis_A81_graphs(k), c) for (k,c) in _tally(p)]...)
      push!(ans, _add_common_root(n, vcat(g...)))
    end
  end
  ans
end

# unlabeled trees
trees(n::Int) = oeis_A55_graphs(n) # alias
function oeis_A55_graphs(n::Int)
  n == 1 && return [SimpleGraph{UInt8}(1)]
  ans = iseven(n) ? [_connect_root(g[1], g[2]) for g in _sym(oeis_A81_graphs(n÷2), 2)] : Graph{Undirected}[]
  for p in _partitions(n-1,n-1,(n-1)÷2)
    for g in Base.Iterators.product([_sym(oeis_A81_graphs(k), c) for (k,c) in _tally(p)]...)
      push!(ans, _add_common_root(n, vcat(g...)))
    end
  end
  ans
end

function _add_common_root(n::Int, gs::AbstractVector{<:Graph})
  g = Graph{Undirected}(n)
  i = 1
  for gi in gs
    add_edge!(g, 1, 1+i)
    for e in edges(gi)
      add_edge!(g, src(e)+i, dst(e)+i)
    end
    i += nv(gi)
  end
  return g
end

function _connect_root(g1::Graph, g2::Graph)
  g = Graph{Undirected}(nv(g1) + nv(g2))
  n1 = nv(g1)
  for e in edges(g1)
    add_edge!(g, src(e), dst(e))
  end
  for e in edges(g2)
    add_edge!(g, src(e)+n1, dst(e)+n1)
  end
  add_edge!(g, 1, 1+n1)
  return g
end

function automorphisms(g::Graph)
  return collect(automorphism_group(g))
  #n = nv(g)
  #ans = Generic.Perm{UInt8}[]
  #p = UInt8.(1:n)
  #for f in automorphism_group(g)
  #  for (i,j) in f
  #    p[i] = j
  #  end
  #  push!(ans, Generic.perm(p[:]))
  #end
  #ans
end

# color a tree with 1:n+1 so that adjacent vertices have different colors
# very naive brute force solution
function colors(g::Graph, n::Int)
  k = nv(g)
  ans = Vector{Int}[]
  for c in Base.Iterators.product(repeat([1:n+1], k)...)
    if all(e -> c[src(e)] != c[dst(e)], edges(g))
      push!(ans, [ci for ci in c])
    end
  end
  ans
end

function multi_trees(n::Int)
  ans = MultiGraph[]
  for m in 2:n+1
    for g in trees(m)
      es = collect(edges(g))
      ess = [sort!([src(e), dst(e)]) for e in es]
      Aut = automorphisms(g)
      edge_configs = _sym(1:m-1, n-m+1)
      while length(edge_configs) > 0
        edge_c = pop!(edge_configs)
        multi = ones(Int, m-1)
        for i in edge_c
          multi[i] += 1
        end
        aut = 0
        for f in Aut
          f_edge_c = sort!([findfirst(==(sort!([f(src(es[i])), f(dst(es[i]))])), ess) for i in edge_c])
          if f_edge_c == edge_c aut += 1 end
          filter!(!=(f_edge_c), edge_configs)
        end
        mul_g = MultiGraph(_copy(g), Dict([e => multi[i] for (i, e) in enumerate(es)]), aut*prod(multi))
        push!(ans, mul_g)
      end
    end
  end
  ans
end

# ###############################################################################
# # 
# # low-level in-place arithmetic operators
# # should probably be added into Nemo upstream
# # 
# function set!(x::fmpz, y::Int)
#   ccall((:fmpz_set_si, Nemo.libflint), Nothing, (Ref{fmpz}, Int), x, y)
#   return x
# end
# function set!(x::QQ, y::QQ)
#   ccall((:fmpz_set, Nemo.libflint), Nothing, (Ref{QQ}, Ref{QQ}), x, y)
#   return x
# end
function Nemo.set!(x::QQFieldElem, n::Int, d::Int)
  ccall((:fmpq_set_si, Nemo.libflint), Nothing, (Ref{QQFieldElem}, Int, UInt), x, n, UInt(d))
  return x
end
# function set!(x::QQ, n::fmpz, d::fmpz)
#   ccall((:QQ_set_fmpz_frac, Nemo.libflint), Nothing, (Ref{QQ}, Ref{fmpz}, Ref{fmpz}), x, n, d)
#   return x
# end
# function Nemo.mul!(z::QQ, x::QQ, y::fmpz)
#   ccall((:QQ_mul_fmpz, Nemo.libflint), Nothing, (Ref{QQ}, Ref{QQ}, Ref{fmpz}), z, x, y)
#   return z
# end
# function Nemo.mul!(z::QQ, x::QQ, y::Int)
#   ccall((:QQ_mul_si, Nemo.libflint), Nothing, (Ref{QQ}, Ref{QQ}, Int), z, x, y)
#   return z
# end
# function Nemo.sub!(z::QQ, x::QQ, y::Int)
#   ccall((:QQ_sub_si, Nemo.libflint), Nothing, (Ref{QQ}, Ref{QQ}, Int), z, x, y)
#   return z
# end
# function Nemo.inv!(x::QQ, y::QQ)
#   ccall((:QQ_inv, Nemo.libflint), Nothing, (Ref{QQ}, Ref{QQ}), x, y)
#   return x
# end
# function pow!(c::QQ, a::QQ, b::Int)
#   iszero(a) && b < 0 && throw(DivideError())
#   ccall((:QQ_pow_si, Nemo.libflint), Nothing, (Ref{QQ}, Ref{QQ}, Int), c, a, b)
#   return c
# end
# function div!(z::QQ, x::QQ, y::QQ)
#   iszero(y) && throw(DivideError())
#   ccall((:QQ_div, Nemo.libflint), Nothing, (Ref{QQ}, Ref{QQ}, Ref{QQ}), z, x, y)
#   return z
# end
# 
# muleq!(x::QQ, y::QQ) = mul!(x, x, y)
# muleq!(x::QQ, y::fmpz) = mul!(x, x, y)
# muleq!(x::QQ, y::Int) = mul!(x, x, y)
# subeq!(x::QQ, y::Int) = Nemo.sub!(x, x, y)
# poweq!(x::QQ, y::Int) = pow!(x, x, y)
# diveq!(x::QQ, y::QQ) = div!(x, x, y)
# Nemo.inv!(x::QQ) = Nemo.inv!(x, x)

###############################################################################
# 
# miscellaneous functions 
# 

function _part(n::T, rem::T, k::T, m::T, part::Vector{T}, ans::Vector{Generic.Partition{T}}) where T <: Integer
  rem == 0 && (push!(ans, Generic.Partition(n, part[1:end-k], false)); return)
  k <= 0 && return
  for v in min(rem, m):-T(1):T(1)
    part[end-k+1] = v
    _part(n, rem-v, k-T(1), v, part, ans)
  end
end
# partitions of n with at most k numbers each <= m
function _partitions(n::T, k::T=n, m::T=n) where T <: Integer
  ans = Generic.Partition{T}[]
  _part(n, n, k, m, Vector{T}(undef, k), ans)
  return ans
end

# similar to combinations but allow multiplicities
function _sym(v::AbstractVector{T}, k::Int) where T
  n = length(v)
  ans = Vector{T}[]
  _sym_dfs!(ans, Vector{T}(undef, k), v, n, k)
  return ans
end
function _sym_dfs!(ans::Vector{Vector{T}}, comb::Vector{T}, v::AbstractVector{T}, n::Int, k::Int) where T
  k < 1 && (pushfirst!(ans, comb[:]); return)
  for m in n:-1:1
    comb[k] = v[m]
    _sym_dfs!(ans, comb, v, m, k - 1)
  end
end

# create a dictionary of counts
function _tally(p::AbstractVector{T}) where T
  ans = Dict{T, Int}()
  for k in p
    ans[k] = k in keys(ans) ? ans[k] + 1 : 1
  end
  ans
end

#=
# for plotting graphs
using GraphPlot
function Base.show(io::IO, mi::MIME"text/html", g::SimpleGraph)
  cm = GraphPlot.Compose.cm
  GraphPlot.Compose.set_default_graphic_size(3cm, 3cm)
  show(io, mi, gplot(g))
end
function Base.show(io::IO, mi::MIME"text/html", g::MultiGraph)
  cm = GraphPlot.Compose.cm
  GraphPlot.Compose.set_default_graphic_size(3cm, 3cm)
  show(io, mi, gplot(g.g, edgelabel=[g[e] > 1 ? g[e] : "" for e in edges(g.g)])) #, nodelabel=collect(1:nv(g.g))))
end
=#

################################################################
### Gromov-Witten Invariants
################################################################

@doc raw"""
     gromov_witten_invariant(d::Int, ns::Vector{Int})
     gromov_witten_invariant(d::Int, ns::Int...)

Return the Gromov-Witten invariant $N_d^{(d_1, \dots, d_k)}$, where $d_1, \dots, d_k$ are the integers given by `ns`.

!!! note
    With respect to notation, we refer to [Dan14](@cite).

# Examples
```jldoctest
julia> [gromov_witten_invariant(d, 5) for d = 1:3]
3-element Vector{QQFieldElem}:
 2875
 4876875//8
 8564575000//27

```

```jldoctest
julia> gromov_witten_invariant(2, 4, 2)
92448

julia> gromov_witten_invariant(2, 3, 3)
423549//8

julia> gromov_witten_invariant(2, 3, 2, 2)
22518

julia> gromov_witten_invariant(2, 2, 2, 2, 2)
9792

```
"""
function gromov_witten_invariant(d::Int, ns::Vector{Int})
  k = length(ns)
  K = kontsevich_moduli_space(k+3, d)
  return integral(K, hypersurface(ns))
end

gromov_witten_invariant(d::Int, ns::Int...) = gromov_witten_invariant(d, collect(ns))


function instanton_numbers(d::Int, ns::Vector{Int})
end

instanton_numbers(d::Int, ns::Int...) = instanton_numbers(d, collect(ns))



