```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Localization and Bott's Formula

Recall that our focus in this chapter is on abstract intersection theory: We discuss computations which manipulate collections of data
referred to as abstract varieties, and we interprete the results as applying to all (smooth projective complex) varieties sharing the data.
The tools presented in this section allow for more efficient computations in the case of varieties with a (split) torus action whose
fixed point set is finite. They are based on localization and a version of Bott's formula which is formulated in the language of
equivariant intersection theory. See [Dan14](@cite) and the references cited there.

Using Bott's formula in enumerative geometry goes back to [ES02](@cite). We quote from that paper:

> Many parameter spaces carry natural actions of algebraic tori, in particular those coming from projective enumerative problems. In 1967, Bott gave a residue formula that allows one to express the degree of certain zero-cycles on a smooth complete variety with an action of an algebraic torus in terms of local contributions supported on the components of the fixpoint set. These components tend to have much simpler structure than the whole space; indeed, in many interesting cases, including all the examples of the present paper, the fixpoints are actually isolated.

We represent an *abstract variety with a torus action* by specifying its dimension together with the fixed points of the action and, possibly, further data.

!!! note
    In order to work with a version of Bott's formula for orbifolds, it is allowed to specify multiplicities at the fixed points. See the section on Kontsevich moduli spaces.

An *abstract equivariant vector bundle under a torus action*  is represented by its rank and its base variety, together with its localizations at the fixed points.

!!! note
    Recall that an equivariant vector bundle over a point is a representation of the group under consideration (in our case, a torus).


## Torus Representations

### Types

For our purposes here, we offer the type `TnRep`.

### Constructors

```@docs
tn_representation(w::Vector{<:IntegerUnion})
```

## Operations on Torus Representations

```@docs
dual(F::TnRep)
```

## Varieties With a Torus Action

### Types

The OSCAR type for abstract varieties with a torus action is `TnVariety`.

### Constructors

```@docs
tn_variety(n::Int, points::Vector{Pair{P, Int}}) where P
```

### Specialized Constructors

```@docs
tn_grassmannian(k::Int, n::Int; weights = :int)
```

```@docs
tn_flag_variety(dims::Int...; weights = :int)
```

### Underlying Data of an Abstract Variety With a Torus Action


```@docs
dim(X::TnVariety)
```

```@docs
fixed_points(X::TnVariety)
```

```@docs
tangent_bundle(X::TnVariety)
```

```@docs
tautological_bundles(X::TnVariety)
```

### Further Data Associated to an Abstract Variety With a Torus Action

As for the type `AbstractVariety`, we have the methods `trivial_line_bundle(X::TnVariety)` (alternatively, `OO(X::TnVariety)`),
`cotangent_bundle(X::TnVariety)`, and `euler_number(X::TnVariety)`. Morever, if `X` is of type `TnVariety`, entering `total_chern_class(X)`
returns the total Chern class of the tangent bundle of `X`. Similarly for entering `chern_class(X, k)`.
	

## Abstract Equivariant Vector Bundles Under a Torus Action

### Types

The OSCAR type for an abstract equivariant vector bundle under a torus action is `TnBundle`.

### Constructors

```@docs
tn_bundle(X::TnVariety, r::Int, f::Function)
```

### Underlying Data of an Equivariant Bundle

If `F` is of type `TnBundle`, then `rank(F)` and `parent(F)` return the rank and the
underlying variety of `F`, respectively. Moreover, we have:

```@docs
localization(F::TnBundle)
```

### Operations on Abstract Equivariant Vector Bundles

```@docs
dual(F::TnBundle)
```

## Chern Classes and Their Integration

In contrast to the varieties of type `AbstractVariety`, there are no associated Chow rings for the varieties of type `TnVariety`.
In order to work with polynomial expressions in the Chern classes of an abstract equivariant vector bundle, Oscar
internally creates an appropriate polynomial ring. We illustrate this in the examples below.

### Types

To work with with polynomial expressions in Chern classes, we offer the type `TnBundleChern`.
 
### Constructors

```@docs
chern_class(F::TnBundle, f::RingElem)
```

```@docs
total_chern_class(F::TnBundle)
```

### Underlying Data of Chern Classes

```@docs
tn_bundle(c::TnBundleChern)
```

```@docs
polynomial(c::TnBundleChern)
```

### Operations on Chern Classes

The usual arithmetic operations are available.

###### Examples

```jldoctest
julia> G = tn_grassmannian(1, 3);

julia> T = tangent_bundle(G)
TnBundle of rank 2 on TnVariety of dim 2 with 3 fixed points

julia> c1 = chern_class(T, 1)
Chern class c[1] of TnBundle of rank 2 on TnVariety of dim 2 with 3 fixed points

julia> c2 = chern_class(T, 2)
Chern class c[2] of TnBundle of rank 2 on TnVariety of dim 2 with 3 fixed points

julia> c = c1^2-3*c2
Chern class c[1]^2 - 3*c[2] of TnBundle of rank 2 on TnVariety of dim 2 with 3 fixed points

julia> typeof(c)
TnBundleChern

```

### Integration

```@docs
integral(c::TnBundleChern)
```

## Examples: Linear Subspaces on Hypersurfaces

```@docs
linear_subspaces_on_hypersurface(k::Int, d::Int; bott::Bool = true)
```

## Kontsevich Moduli Spaces

## Examples: Rational Curves on Complete Intersections

