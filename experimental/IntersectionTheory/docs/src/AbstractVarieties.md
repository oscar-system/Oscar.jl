```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Abstract varieties

An *abstract variety* $X$ of dimension $n$ is determined by its numerical Chow ring
$\mathrm{N}^*(X)_{\mathbb Q}$, together with an integration map (degree map)
$\int_X\colon \mathrm{N}^n(X)_{\mathbb Q}\to \mathbb Q$. It may additionally carry
a tangent bundle, a polarization, tautological bundles, and a structure map to a base variety.

Abstract varieties can be constructed either from scratch by specifying a graded ring
and integration map, or via specialized constructors for standard algebraic-geometric
objects such as projective spaces, Grassmannians, flag varieties, complete intersections,
and projective bundles.

## Types

The OSCAR type for abstract varieties is `AbstractVariety`.

## Constructors

```@docs
abstract_variety(n::Int, A::MPolyDecRingOrQuo)
```

### Specialized constructors


```@docs
abstract_point(; base::Ring=QQ)
```

```@docs
abstract_projective_space(n::Int; base::Ring = QQ, symbol::String = "h")
```

```@docs
abstract_grassmannian(k::Int, n::Int; bott::Bool = false, weights = :int, base::Ring = QQ, symbol::String = "c")
```

```@docs
abstract_flag_variety(dims::Int...; base::Ring = QQ, symbol::String = "c")
```

### New varieties from given varieties/bundles

```@docs
complete_intersection(X::AbstractVariety, degs::Int...)
```

```@docs
projective_bundle(F::AbstractBundle; symbol::String = "z")
```

```@docs
abstract_hirzebruch_surface(n::Int)
```

```@docs
flag_bundle(F::AbstractBundle, dims::Int...; symbol::String = "c")
```

```@docs
zero_locus_section(F::AbstractBundle; class::Bool = false)
```

```@docs
degeneracy_locus(F::AbstractBundle, G::AbstractBundle, k::Int; class::Bool=false)
```

!!! note
    Products and blow-ups are described elsewhere.

## Underlying data of an abstract variety

An abstract variety is made up from (a selection of) the data discussed here:

```@docs
dim(X::AbstractVariety)
```

```@docs
chow_ring(X::AbstractVariety)
```

```@docs
base(X::AbstractVariety)
```

```@docs
point_class(X::AbstractVariety)
```

```@docs
tangent_bundle(X::AbstractVariety)
```

```@docs
polarization(X::AbstractVariety)
```

```@docs
tautological_bundles(X::AbstractVariety)
```

```@docs
structure_map(X::AbstractVariety)
```

## Further data associated to an abstract variety


```@docs
trivial_line_bundle(X::AbstractVariety)
```

```@docs
line_bundle(X::AbstractVariety, n::RingElement)
```

```@docs
cotangent_bundle(X::AbstractVariety)
```

```@docs
canonical_class(X::AbstractVariety)
```

```@docs
canonical_bundle(X::AbstractVariety)
```

```@docs
degree(X::AbstractVariety)
```

```@docs
hilbert_polynomial(X::AbstractVariety)
```

```@docs
basis(X::AbstractVariety)
```

```@docs
betti_numbers(X::AbstractVariety)
```

```@docs
intersection_matrix(X::AbstractVariety)
```

```@docs
dual_basis(X::AbstractVariety)
```

```@docs
euler_number(X::AbstractVariety)
```

!!! note
    If `X` is of type `AbstractVariety`, entering `total_chern_class(X)` returns the total Chern class of the tangent bundle of `X`.
    Similarly for entering `chern_class(X, k)`, `todd_class(X)`, `total_pontryagin_class(X)`, and `pontryagin_class(X, k)`.
    Moreover, `gens(X)` returns the generators of the Chow ring of `X`.

## Operations on abstract varieties

```@docs
product(X::AbstractVariety, Y::AbstractVariety)
```

!!! note
    Blow-Ups are described in their own section.

## Integrating Chow ring elements

```@julia
integral(c::Union{MPolyDecRingElem, MPolyQuoRingElem})
```

Given an element `c` of the Chow ring of an abstract variety, return the integral of `c`.

!!! note
    If the abstract variety has been given a (unique) point class,
    then the integral will be an element of the coefficient ring of the Chow ring.
    That is, typically, in the applications we discuss here, it will be a rational number (the degree of the 0-dimensional part
    of `c`) or an element of a function field of type $\mathbb Q(t_1, \dots, t_r)$. If one of the conditions is not fulfilled, the 0-dimensional
    part of `c` is returned.


###### Examples

```jldoctest
julia> G = abstract_grassmannian(2, 5)
AbstractVariety of dim 6

julia> p = point_class(G)
c[2]^3

julia> integral(p)
1

```

```jldoctest
julia> T, (t, ) = polynomial_ring(QQ, [:t])
(Multivariate polynomial ring in 1 variable over QQ, QQMPolyRingElem[t])

julia> QT = fraction_field(T)
Fraction field
  of multivariate polynomial ring in 1 variable over QQ

julia> P3 = abstract_projective_space(3, base = QT)
AbstractVariety of dim 3

julia> h = gens(P3)[1]
h

julia> integral(t^2*h^3+t*h)
t^2

```

*Lines on a general cubic hypersurface in P^3:*

```jldoctest
julia> G = abstract_grassmannian(2, 4)
AbstractVariety of dim 4

julia> Q = tautological_bundles(G)[2]
AbstractBundle of rank 2 on AbstractVariety of dim 4

julia> E = symmetric_power(Q, 3)
AbstractBundle of rank 4 on AbstractVariety of dim 4

julia> integral(top_chern_class(E))
27

```

*Lines on a general complete intersection Calabi-Yau threefold of type (2,2,2,2):*

```jldoctest
julia> G = abstract_grassmannian(2, 4+4)
AbstractVariety of dim 12

julia> S = tautological_bundles(G)[1]
AbstractBundle of rank 2 on AbstractVariety of dim 12

julia> E = symmetric_power(S, 2)
AbstractBundle of rank 3 on AbstractVariety of dim 12

julia> integral(top_chern_class(E)^4)
512

```
