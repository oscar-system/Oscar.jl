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

@doc raw"""
    dual(F::TnRep)
    det(F::TnRep)
    +(F::TnRep, G::TnRep)
    *(F::TnRep, G::TnRep)

Return the dual of `F`, the determinant of `F`, the sum `F` $+$ `G`, and the tensor product of `F` and `G`, respectively.

# Examples
```jldoctest
julia> F = tn_representation([1, 2, 3])
TnRep(3, ZZRingElem[1, 2, 3])

julia> G = tn_representation([1, 1, 1])
TnRep(3, ZZRingElem[1, 1, 1])

julia> F*G
TnRep(9, ZZRingElem[2, 2, 2, 3, 3, 3, 4, 4, 4])

```
"""
dual(F::TnRep) = TnRep(-F.w)
det(F::TnRep) = TnRep([sum(F.w)])
+(F::TnRep, G::TnRep) = TnRep(vcat(F.w, G.w))
*(F::TnRep, G::TnRep) = TnRep([a+b for a in F.w for b in G.w])

top_chern_class(F::TnRep) = prod(F.w)
function chern_class(F::TnRep, n::Int)
  sum(prod(F.w[i] for i in c) for c in combinations(F.n, n))
end
function _sym(k::Int, n::Int)
  k == 0 && return [Int[]]
  vcat([[push!(c, i) for c in _sym(k-1,i)] for i in 1:n]...)
end

@doc raw"""
    tn_representation(w::Vector{<:IntegerUnion})

Return the representation of a (split) torus of rank `length(w)` specified by the weights `w`.

# Examples
```jldoctest
julia> tn_representation([1, 2, 3])
TnRep(3, ZZRingElem[1, 2, 3])

```
"""
tn_representation(w::Vector{<:IntegerUnion}) = TnRep(w)


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
    localization(F::TnBundle)

Return the localization of `F` at the fixed points of the given torus action.

# Examples
```jldoctest
julia> G = tn_grassmannian(1, 3);

julia> T = tangent_bundle(G)
TnBundle of rank 2 on TnVariety of dim 2 with 3 fixed points

julia> V = fixed_points(G)
3-element Vector{Pair{Vector{Int64}, Int64}}:
 [1] => 1
 [2] => 1
 [3] => 1

julia> f = localization(T);

julia> P1 = V[1][1];

julia> f(P1)
TnRep(2, ZZRingElem[1, 2])

julia> P2 = V[2][1];

julia> f(P2)
TnRep(2, ZZRingElem[-1, 1])

julia> P3 = V[3][1];

julia> f(P3)
TnRep(2, ZZRingElem[-2, -1])

```
"""
localization(F::TnBundle) = F.loc

@doc raw"""
    TnVariety(n::Int, points)

The type of an abstract variety with a (split) torus action, represented by its dimension `n` and a `Vector` determining the fixed points of the action together with their multiplicities.
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

Base.show(io::IO, X::TnVariety) = print(io,
  "TnVariety of dim ", X.dim, " with ", length(X.points), " fixed points")

@doc raw"""
    tn_bundle(X::TnVariety, r::Int, f::Function)

Return an abstract equivariant vector bundle on `X` by specifying its rank together with a function which gives the localization of the vector bundle at each fixed point of the torus action on `X`.
"""
tn_bundle(X::TnVariety, r::Int, f::Function) = TnBundle(X, r, f) 

@doc raw"""
    tn_variety(n::Int, points::Vector{Pair{P, Int}}) where P

Return an abstract variety with a (split) torus action, represented by its dimension `n` and a `Vector` specifying the fixed points of the action together with their multiplicities.

!!! note
    Specifying multiplicities at the fixed points allows one to work with a version of Bott's formula for orbifolds. See the section on Kontsevich moduli spaces in the documentation.
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

Return the fixed points representing `X` and their multiplicities.

!!! note
    Specifying multiplicities at the fixed points allows one to work with a version of Bott's formula for orbifolds. See the section of the documentation on Kontsevich moduli spaces.

# Examples
```jldoctest
julia> G = tn_grassmannian(2, 5);

julia> V = fixed_points(G)
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

julia> P = V[10][1]
2-element Vector{Int64}:
 4
 5

```
"""
fixed_points(X::TnVariety) = X.points

@doc raw"""
     tangent_bundle(X::TnVariety)

Return the tangent bundle of `X`.

# Examples
```jldoctest
julia> G = tn_grassmannian(1, 3);

julia> T = tangent_bundle(G)
TnBundle of rank 2 on TnVariety of dim 2 with 3 fixed points

julia> V = fixed_points(G)
3-element Vector{Pair{Vector{Int64}, Int64}}:
 [1] => 1
 [2] => 1
 [3] => 1

julia> f = localization(T);

julia> P1 = V[1][1];

julia> f(P1)
TnRep(2, ZZRingElem[1, 2])

julia> P2 = V[2][1];

julia> f(P2)
TnRep(2, ZZRingElem[-1, 1])

julia> P3 = V[3][1];

julia> f(P3)
TnRep(2, ZZRingElem[-2, -1])

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
 TnBundle of rank 2 on TnVariety of dim 6 with 10 fixed points
 TnBundle of rank 3 on TnVariety of dim 6 with 10 fixed points

```
"""
tautological_bundles(X::TnVariety) = X.bundles  ### CHECK bundles( --> tautological_bundles(

euler_number(X::TnVariety) = sum(1//ZZ(e) for (p,e) in X.points) # special case of Bott's formula
cotangent_bundle(X::TnVariety) = dual(X.T)
trivial_line_bundle(X::TnVariety) = TnBundle(X, 1, p -> TnRep([0]))
(OO)(X::TnVariety) = trivial_line_bundle(X)

@doc raw"""
    dual(F::TnBundle)
    symmetric_power(F::TnBundle, k::Int)
    exterior_power(F::TnBundle, k::Int)
    det(F::TnBundle)
    +(F::TnBundle, G::TnBundle)
    *(F::TnBundle, G::TnBundle)

Return the dual of `F`, the `k`-th symmetric power of `F`, the `k`-th exterior power of `F`, the determinant of `F`, the sum `F` $+$ `G`, and the tensor product of `F` and `G`, respectively.
"""
dual(F::TnBundle) = TnBundle(F.parent, F.rank, p -> dual(F.loc(p)))
det(F::TnBundle) = TnBundle(F.parent, 1, p -> det(F.loc(p)))
+(F::TnBundle, G::TnBundle) = TnBundle(F.parent, F.rank + G.rank, p -> F.loc(p) + G.loc(p))
*(F::TnBundle, G::TnBundle) = TnBundle(F.parent, F.rank * G.rank, p -> F.loc(p) * G.loc(p))

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

@doc raw"""
     total_chern_class(F::TnBundle)
     chern_class(F::TnBundle, k::Int)
     top_chern_class(F::TnBundle)

Return the total Chern class, the `k`-th Chern class, and the top Chern class of `F`, respectively.

# Examples

```jldoctest
julia> G = tn_grassmannian(1, 3);

julia> T = tangent_bundle(G)
TnBundle of rank 2 on TnVariety of dim 2 with 3 fixed points

julia> total_chern_class(T::TnBundle)
Chern class c[1] + c[2] + 1 of TnBundle of rank 2 on TnVariety of dim 2 with 3 fixed points

```
"""
total_chern_class(F::TnBundle) = TnBundleChern(F, 1+sum(gens(_get_ring(F))))
chern_class(F::TnBundle, k::Int) = TnBundleChern(F, total_chern_class(F).c[k])
top_chern_class(F::TnBundle) = chern_class(F, F.rank)

@doc raw"""
     chern_class(F::TnBundle, f::RingElem)

Return the evaluation of `f` in the Chern classes of `F`.

# Examples

```jldoctest
julia> G = tn_grassmannian(1, 3);

julia> T = tangent_bundle(G)
TnBundle of rank 2 on TnVariety of dim 2 with 3 fixed points

julia> R, (x, y) = polynomial_ring(QQ, [:x, :y])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> f = x*y^2
x*y^2

julia> c = chern_class(T, f)
Chern class c[1]*c[2]^2 of TnBundle of rank 2 on TnVariety of dim 2 with 3 fixed points

julia> parent(polynomial(c))
Multivariate polynomial ring in 2 variables over QQ graded by
  c[1] -> [1]
  c[2] -> [2]

```
"""
chern_class(F::TnBundle, x::RingElem) = begin
  R = _get_ring(F)
  @assert length(gens(R)) == length(gens(parent(x)))
  TnBundleChern(F, evaluate(x, gens(R)))
end

@doc raw"""
     tn_bundle(c::TnBundleChern)

Return the `tn_bundle` to which `c` belongs.
"""
tn_bundle(c::TnBundleChern) = c.F

@doc raw"""
    polynomial(c::TnBundleChern)

Return the polynomial representing `c`.

# Examples
```jldoctest
julia> G = tn_grassmannian(1, 3);

julia> T = tangent_bundle(G)
TnBundle of rank 2 on TnVariety of dim 2 with 3 fixed points

julia> f = polynomial(total_chern_class(T))
c[1] + c[2] + 1

julia> parent(f)
Multivariate polynomial ring in 2 variables over QQ graded by
  c[1] -> [1]
  c[2] -> [2]

```
"""
polynomial(c::TnBundleChern) = c.c

total_chern_class(X::TnVariety) = total_chern_class(X.T)
chern_class(X::TnVariety, k::Int) = chern_class(X.T, k)

@doc raw"""
     integral(c::TnBundleChern)

Return the integral of `c`.

# Examples
```jldoctest
julia> G = tn_grassmannian(2, 4)
TnVariety of dim 4 with 6 fixed points

julia> Q = tautological_bundles(G)[2];

julia> E = symmetric_power(Q, 3)
TnBundle of rank 4 on TnVariety of dim 4 with 6 fixed points

julia> integral(top_chern_class(E))
27

```
"""
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

Return the Grassmannian $\mathrm{G}(k, n)$ of `k`-dimensional subspaces of an 
`n`-dimensional standard vector space as a `TnVariety`, where the action is induced by
the diagonal action with `weights` on the standard vector space.

!!! note
    The fixed points correspond to the ${n}\choose{k}$ coordinate subspaces of dimension $k$ in the standard vector space.

# Examples
```jldoctest
julia> G = tn_grassmannian(3,5)  # all weights are 1
TnVariety of dim 6 with 10 fixed points

julia> V = fixed_points(G)
10-element Vector{Pair{Vector{Int64}, Int64}}:
 [1, 2, 3] => 1
 [1, 2, 4] => 1
 [1, 3, 4] => 1
 [2, 3, 4] => 1
 [1, 2, 5] => 1
 [1, 3, 5] => 1
 [2, 3, 5] => 1
 [1, 4, 5] => 1
 [2, 4, 5] => 1
 [3, 4, 5] => 1

julia> P = V[10][1]
3-element Vector{Int64}:
 3
 4
 5

```
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

Given integers, say, $d_1, \dots, d_{k}, n$ with $0 < d_1 < \dots < d_{k} < n$, 
return the abstract flag variety $\mathrm{F}(d_1, \dots, d_{k}; n)$ of nested sequences of subspaces of 
dimensions $d_1, \dots, d_{k}$ of an $n$-dimensional standard vector space as a `TnVariety`, where the 
action is induced by the diagonal action with `weights` on the standard vector space.

# Examples
```jldoctest
julia> F = tn_flag_variety(1,3,4)  # all weights are 1 
TnVariety of dim 5 with 12 fixed points

julia> fixed_points(F)
12-element Vector{Pair{Vector{Vector{Int64}}, Int64}}:
 [[1], [2, 3], [4]] => 1
 [[1], [2, 4], [3]] => 1
 [[1], [3, 4], [2]] => 1
 [[2], [1, 3], [4]] => 1
 [[2], [1, 4], [3]] => 1
 [[2], [3, 4], [1]] => 1
 [[3], [1, 2], [4]] => 1
 [[3], [1, 4], [2]] => 1
 [[3], [2, 4], [1]] => 1
 [[4], [1, 2], [3]] => 1
 [[4], [1, 3], [2]] => 1
 [[4], [2, 3], [1]] => 1

```
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
