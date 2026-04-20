```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Abstract varieties

An *abstract variety* $X$ of dimension $n$ is determined by its (numerical) Chow ring
$\mathrm{N}^*(X)_{\mathbb Q} = \bigoplus^n_{c=0}\mathrm{N}^c(X)_{\mathbb Q}$,
together with an integration map (degree map) $\int_X\colon \mathrm{N}^n(X)_{\mathbb Q}\to \mathbb Q$.

Each Chow ring $\mathrm{N}^*(X)_{\mathbb Q}$ is given by finitely many generators and relations.
That is, it is implemented as the quotient of a $\mathbb Z$-graded multivariate polynomial ring
over $\mathbb Q$ modulo a homogeneous ideal. It has Krull dimension zero and is, thus, a
finite-dimensional $\mathbb Q$-vector space. The Betti numbers
$\beta_c(X) = \dim_{\mathbb Q} \mathrm{N}^c(X)_{\mathbb Q}$
satisfy the relations $\beta_c(X) = \beta_{n-c}(X)$ for each $c$. In particular,
$\beta_n = \beta_0 = 1$.

To specify an integration map means to specify a point class, that is, a (unique) degree-$n$ element
of the Chow ring that integrates to one. Additionally, an abstract variety may carry a tangent bundle, a
polarization, tautological bundles, and a structure map to a base variety. See the setter
functions in section [Some particular constructions](@ref).

Abstract varieties can be constructed either from scratch by specifying a graded ring
as above and a point class, or via specialized constructors for standard algebraic-geometric
objects such as projective spaces, Grassmannians, flag varieties, complete intersections,
and projective bundles.

!!! note
    The constructor `abstract_variety` discussed below gives the expert user some freedom when
    constructing an object of type `AbstractVariety`. It allows one, for example, to start from the
    underlying graded polynomial ring of the Chow ring, and add its defining relations step by step.
    In fact, not all applications require that we specify all relations. See section [Some Particular Constructions](@ref)
	for an example where the top-dimensional part of the constructed ring is more than 1-dimensional.

## Types

The OSCAR type for abstract varieties is `AbstractVariety`.

## Constructors

### General varieties

```@docs
abstract_variety(n::Int, A::MPolyDecRingOrQuo)
```

```@docs
abstract_point(; base::Ring=QQ)
```

```@docs
abstract_curve(g::IntegerUnion; base::Ring = QQ)
```

### Flag bundles

The concept of flag bundles provides fundamental classifying spaces in enumerative geometry.
In Oscar, abstract flag bundles are constructed using the function `flag_bundle` which, in particular,
allows one to implement abstract projective spaces, Grassmannians, flag varieties, and projective bundles.
In addition, there are specialized constructors for the latter varieties some of which which rely on
different recipes for representing Chow rings  in terms of generators and relations.

```@docs
abstract_projective_space(n::Int; base::Ring = QQ, symbol::String = "h")
```

```@docs
abstract_grassmannian(k::Int, n::Int; bott::Bool = false, weights = :int, base::Ring = QQ, symbol::String = "c")
```

```@docs
abstract_flag_variety(dims::Int...; base::Ring = QQ, symbol::String = "c")
```

```@docs
projective_bundle(F::AbstractBundle; symbol::String = "z")
```

```@docs
flag_bundle(F::AbstractBundle, dims::Int...; symbol::String = "c")
```

### Further standard varieties

```@docs
complete_intersection(X::AbstractVariety, degs::Int...)
```

```@docs
abstract_quadric(n::Int; base::Ring = QQ)
```

```@docs
zero_locus_section(F::AbstractBundle; class::Bool = false)
```

```@docs
degeneracy_locus(F::AbstractBundle, G::AbstractBundle, k::Int; class::Bool=false)
```

### Special varieties

```@docs
abstract_K3_surface(g::IntegerUnion; base::Ring = QQ)
```

```@docs
abstract_hirzebruch_surface(n::Int)
```

```@docs
abstract_cayley_plane(; base::Ring = QQ)
```

```@docs
abstract_cayley_grassmannian(; base::Ring = QQ)
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
    [Blow-Ups](@ref) are described in their own section.

## Integrating Chow ring elements

```@julia
integral(c::Union{MPolyDecRingElem, MPolyQuoRingElem})
```

Given an element `c` of the Chow ring of an abstract variety, say, `X`, return the integral of `c`.

!!! note
    If `X` has been given a point class, and the top Betti number of  `X` is 1,
    then the integral will be an element of the coefficient ring of the Chow ring.
    That is, in the applications we discuss here, it will be a rational number (the degree of the 0-dimensional part
    of `c`) or an element of a function field of type $\mathbb Q(t_1, \dots, t_r)$. If not both conditions on `X` are
    fulfilled, the 0-dimensional part of `c` will be returned.


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
