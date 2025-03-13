###############################################################################
#
# TnRep - n-dim representation of a torus, specified by its weights
#
@doc raw"""
    TnRep(w::Vector)

The type of a representation of a torus, specified by its weights.
"""
struct TnRep
  n::Int
  w::Vector
  function TnRep(w::Vector{W}) where W
    # be sure to use ZZRingElem to avoid overflow
    W == Int && return new(length(w), ZZ.(w))
    new(length(w), w)
  end
end
dual(F::TnRep) = TnRep(-F.w)
+(F::TnRep, G::TnRep) = TnRep(vcat(F.w, G.w))
*(F::TnRep, G::TnRep) = TnRep([a+b for a in F.w for b in G.w])
det(F::TnRep) = TnRep([sum(F.w)])
top_chern_class(F::TnRep) = prod(F.w)
function chern_class(F::TnRep, n::Int)
  sum(prod(F.w[i] for i in c) for c in combinations(F.n, n))
end
function _sym(k::Int, n::Int)
  k == 0 && return [Int[]]
  vcat([[push!(c, i) for c in _sym(k-1,i)] for i in 1:n]...)
end

###############################################################################
#
# TnBundle, TnVariety - varieties with a torus action and equivariant bundles
#
# A Tⁿ-abstract_variety X is represented as the set of fixed points X.points, each
# labeled using some value of type P (e.g. an array), and has a multiplicity e
# (orbifold multiplicity);
# 
# A Tⁿ-equivariant bundle on X is represented by its localization/restriction
# to each of the points in X.points, which will be of type `TnRep`.
# They are stored as a function to allow lazy evaluation: this is crucial for
# large examples, since otherwise we may run into memory problems.
abstract type TnVarietyT{P} <: Variety end

@doc raw"""
    TnBundle(X::TnVariety, r::Int, f::Function)

The type of a torus-equivariant bundle, represented by its localizations to the
fixed points of the base abstract_variety.
"""
@attributes mutable struct TnBundle{P, V <: TnVarietyT{P}} <: Bundle
  parent::V
  rank::Int
  loc::Function

  function TnBundle(X::V, r::Int) where V <: TnVarietyT
    P = V.parameters[1]
    new{P, V}(X, r)
  end
  function TnBundle(X::V, r::Int, f::Function) where V <: TnVarietyT
    P = V.parameters[1]
    new{P, V}(X, r, f)
  end
end

@doc raw"""
    TnVariety(n::Int, points)

The type of a abstract_variety with a torus action, represented by the fixed points.
"""
@attributes mutable struct TnVariety{P} <: TnVarietyT{P}
  dim::Int
  points::Vector{Pair{P, Int}}
  T::TnBundle
  bundles::Vector{TnBundle}

  function TnVariety(n::Int, points::Vector{Pair{P, Int}}) where P
    new{P}(n, points)
  end
end

@doc raw"""
    tn_variety(n::Int, points::Vector{Pair{P, Int}}) where P

Return an abstract_variety with a torus action, represented by the fixed points.
"""
tn_variety(n::Int, points::Vector{Pair{P, Int}}) where P = TnVariety(n::Int, points::Vector{Pair{P, Int}})

@doc raw"""
     dim(X::TnVariety)

Return the dimension of `X`.

# Examples
```jldoctest
julia> G = tn_grassmannian(2, 5);

julia> dim(G)
6

```
"""
dim(X::TnVariety) = X.dim

@doc raw"""
     fixed_points(X::TnVariety)

Return the fixed points representing `X`.

# Examples
```jldoctest
julia> G = tn_grassmannian(2, 5);

julia> fixed_points(G)
10-element Vector{Pair{Vector{Int64}, Int64}}:
 [1, 2] => 1
 [1, 3] => 1
 [2, 3] => 1
 [1, 4] => 1
 [2, 4] => 1
 [3, 4] => 1
 [1, 5] => 1
 [2, 5] => 1
 [3, 5] => 1
 [4, 5] => 1

```
"""
fixed_points(X::TnVariety) = X.points

@doc raw"""
     tangent_bundle(X::TnVariety)

Return the tangent bundle of `X`.

# Examples
```jldoctest
julia> G = tn_grassmannian(2, 5);

julia> tangent_bundle(G)
TnBundle of rank 6 on TnVariety of dim 6

```
"""
tangent_bundle(X::TnVariety) = X.T

@doc raw"""
     tautological_bundles(X::TnVariety)

If `X` has been given tautological bundles, return these bundles.

# Examples
```jldoctest
julia> G = tn_grassmannian(2, 5);

julia> tautological_bundles(G)
2-element Vector{TnBundle}:
 TnBundle of rank 2 on TnVariety of dim 6
 TnBundle of rank 3 on TnVariety of dim 6

```
"""
tautological_bundles(X::TnVariety) = X.bundles  ### CHECK bundles( --> tautological_bundles(

euler(X::TnVariety) = sum(1//ZZ(e) for (p,e) in X.points) # special case of Bott's formula
cotangent_bundle(X::TnVariety) = dual(X.T)
trivial_line_bundle(X::TnVariety) = TnBundle(X, 1, p -> TnRep([0]))

dual(F::TnBundle) = TnBundle(F.parent, F.rank, p -> dual(F.loc(p)))
+(F::TnBundle, G::TnBundle) = TnBundle(F.parent, F.rank + G.rank, p -> F.loc(p) + G.loc(p))
*(F::TnBundle, G::TnBundle) = TnBundle(F.parent, F.rank * G.rank, p -> F.loc(p) * G.loc(p))
det(F::TnBundle) = TnBundle(F.parent, 1, p -> det(F.loc(p)))

# avoid computing `_sym` for each F.loc(p)
function symmetric_power(F::TnBundle, k::Int)
  l = _sym(k, F.rank)
  TnBundle(F.parent, binomial(F.rank+k-1, k), p -> (
    Fp = F.loc(p);
    TnRep([sum(Fp.w[i] for i in c) for c in l])))
end
function exterior_power(F::TnBundle, k::Int)
  l = combinations(F.rank, k)
  TnBundle(F.parent, binomial(F.rank, k), p -> (
    Fp = F.loc(p);
    TnRep([sum(Fp.w[i] for i in c) for c in l])))
end

# we want the same syntax `integral(total_chern_class(F))` as in Schubert calculus
# the following ad hoc type represents a formal expression in chern classes of a bundle F
struct TnBundleChern
  F::TnBundle
  c::MPolyDecRingElem
end
for O in [:(+), :(-), :(*)]
  @eval $O(a::TnBundleChern, b::TnBundleChern) = (
    @assert a.F == b.F;
    TnBundleChern(a.F, $O(a.c, b.c)))
end
^(a::TnBundleChern, n::Int) = TnBundleChern(a.F, a.c^n)
*(a::TnBundleChern, n::RingElement) = TnBundleChern(a.F, a.c*n)
*(n::RingElement, a::TnBundleChern) = TnBundleChern(a.F, n*a.c)
Base.show(io::IO, c::TnBundleChern) = print(io, "Chern class $(c.c) of $(c.F)")

# create a ring to hold the chern classes of F
function _get_ring(F::TnBundle)
  if get_attribute(F, :R) === nothing
    r = min(F.parent.dim, F.rank)
    R, _ = graded_polynomial_ring(QQ, :c => 1:r; weights = 1:r)
    set_attribute!(R, :abstract_variety_dim => F.parent.dim)
    set_attribute!(F, :R => R)
  end
  get_attribute(F, :R)
end

total_chern_class(F::TnBundle) = TnBundleChern(F, 1+sum(gens(_get_ring(F))))
total_chern_class(X::TnVariety) = total_chern_class(X.T)
chern_class(F::TnBundle, k::Int) = TnBundleChern(F, total_chern_class(F).c[k])
chern_class(X::TnVariety, k::Int) = chern_class(X.T, k)
top_chern_class(F::TnBundle) = chern_class(F, F.rank)
chern_class(F::TnBundle, x::RingElem) = begin
  R = _get_ring(F)
  @assert length(gens(R)) == length(gens(parent(x)))
  TnBundleChern(F, evaluate(x, gens(R)))
end

function integral(c::TnBundleChern)
  F, R = c.F, parent(c.c)
  X = F.parent
  n, r = X.dim, length(gens(R))
  top = c.c[n].f
  top == 0 && return QQ()
  exp_vec = sum(AbstractAlgebra.exponent_vectors(top))
  idx = filter(i -> exp_vec[i] > 0, 1:r)
  ans = 0
  for (p,e) in X.points # Bott's formula
    Fp = F.loc(p)
    # this avoids the computations of Chern classes that are not needed
    cherns = [i in idx ? chern_class(Fp, i) : QQ() for i in 1:r]
    ans += top(cherns...) * (1 // (e * top_chern_class(X.T.loc(p))))
  end
  ans
end

###############################################################################
#
# Grassmannians and flag varieties
#
# utility function that parses the weight specification

function _parse_weight(n::Int, w)
  w == :int && return ZZ.(collect(1:n))
  w == :poly && return polynomial_ring(QQ, "u#" => 1:n)[2]
  if (w isa AbstractUnitRange) w = collect(w) end
  w isa Vector && length(w) == n && return w
  error("incorrect specification for weights")
end

@doc raw"""
    tn_grassmannian(k::Int, n::Int; weights = :int)

Return a Grassmannian with a torus action, ... .
"""
function tn_grassmannian(k::Int, n::Int; weights = :int)
  @assert k < n
  points = [p=>1 for p in combinations(n, k)]
  d = k*(n-k)
  G = TnVariety(d, points)
  w = _parse_weight(n, weights)
  S = TnBundle(G, k, p -> TnRep([w[i] for i in p]))
  Q = TnBundle(G, n-k, p -> TnRep([w[i] for i in setdiff(1:n, p)]))
  G.bundles = [S, Q]
  G.T = dual(S) * Q
  set_attribute!(G, :description => "Grassmannian Gr($k, $n)")
  return G
end

@doc raw"""
    tn_flag_variety(dims::Int...; weights = :int)

Return a flag variety with a torus action, ... .
"""
function tn_flag_variety(dims::Int...; weights = :int)
  return tn_flag_variety(collect(dims), weights  = weights)
end

function tn_flag_variety(dims::Vector{Int}; weights = :int)
  n, l = dims[end], length(dims)
  ranks = pushfirst!([dims[i+1]-dims[i] for i in 1:l-1], dims[1])
  @assert all(>(0), ranks)
  d = sum(ranks[i] * sum(dims[end]-dims[i]) for i in 1:l-1)
  function enum(i::Int, rest::Vector{Int})
    i == l && return [[rest]]
    [pushfirst!(y, x) for x in combinations(rest, ranks[i]) for y in enum(i+1, setdiff(rest, x))]
  end
  points = [p=>1 for p in enum(1, collect(1:n))]
  Fl = TnVariety(d, points)
  w = _parse_weight(n, weights)
  Fl.bundles = [TnBundle(Fl, r, p -> TnRep([w[j] for j in p[i]])) for (i, r) in enumerate(ranks)]
  Fl.T = sum(dual(Fl.bundles[i]) * sum([Fl.bundles[j] for j in i+1:l]) for i in 1:l-1)
  set_attribute!(Fl, :description => "Flag abstract_variety Flag$(tuple(dims...))")
  return Fl
end

@doc raw"""
    linear_subspaces_on_hypersurface(k::Int, d::Int; bott::Bool = true)

If $n=\frac1{k+1}\binom{d+k}d+k$ is an integer, return the number of $k$-dimensional subspaces 
on a generic hypersurface of degree $d$ in a projective space of dimension $n$.

    lines_on_hypersurface(n::Int; bott::Bool = true)

Return the number of lines on a hypersurface of degree $d = 2*n-3$ in a projective space of dimension $n$.

!!! note
    The function relies on Bott's formula by default. Use `bott = false` to switch to Schubert calculus.

# Examples

```jldoctest
julia> linear_subspaces_on_hypersurface(1,3)
27

julia> linear_subspaces_on_hypersurface(2,5)
420760566875

```

```jldoctest
julia> [lines_on_hypersurface(n) for n=2:10]
9-element Vector{QQFieldElem}:
 1
 27
 2875
 698005
 305093061
 210480374951
 210776836330775
 289139638632755625
 520764738758073845321

```
"""
function linear_subspaces_on_hypersurface(k::Int, d::Int; bott::Bool = true)
  n = binomial(d+k, d) // (k+1)
  is_integer(n) || error("binomial(d+k, d) // (k+1) is not an integer")
  n = Int(n)+k
  if bott == true
    G = tn_grassmannian(k+1, n+1)
  else
    G = abstract_grassmannian(k+1, n+1)
  end
  S, Q = G.bundles
  integral(top_chern_class(symmetric_power(dual(S), d)))
end

function lines_on_hypersurface(n::Int; bott::Bool = true)
  n >= 2 || error("n must be at least 2")
  return linear_subspaces_on_hypersurface(1, 2*n-3)
end
