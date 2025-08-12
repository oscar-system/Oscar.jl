```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```


# Projective schemes

Let ``A`` be a commutative noetherian base ring and
``S = A[x_0,\dots, x_n]`` the standard graded polynomial ring
over ``A``. Then ``X = \mathrm{Proj}(S) = \mathbb P^n_A`` is a
(relative) projective scheme over ``\mathrm{Spec}(A)``.
Similarly, for a homogeneous ideal ``I \subset S`` we have
``X = \mathrm{Proj}(S/I) \subset \mathbb P^n_A`` a (relative)
projective scheme over ``\mathrm{Spec}(A)`` which is a closed
subscheme of ``\mathbb P^n_A`` in a natural way. The majority
of applications will be in the setting where ``A = \mathbb k`` is a
field, but be aware that we also support different base rings
such as the usual four `MPolyRing`, `MPolyQuoRing`, `MPolyLocRing`,
and `MPolyQuoLocRing`.

## Abstract types and basic interface
The abstract type for such projective schemes is
```julia
AbsProjectiveScheme{CoeffRingType, RingType} where {CoeffRingType<:Ring}
```
where, in the above notation, `CoeffRingType` denotes the type of `A`
and `RingType` the type of either `S` or `S/I`, respectively.
The abstract type comes with the following interface:
```@docs
base_ring(X::AbsProjectiveScheme)
base_scheme(X::AbsProjectiveScheme)
homogeneous_coordinate_ring(P::AbsProjectiveScheme)
relative_ambient_dimension(X::AbsProjectiveScheme)
ambient_coordinate_ring(P::AbsProjectiveScheme)
ambient_space(P::AbsProjectiveScheme)
defining_ideal(X::AbsProjectiveScheme)
affine_cone(X::AbsProjectiveScheme)
homogeneous_coordinates_on_affine_cone(X::AbsProjectiveScheme)
covered_scheme(P::AbsProjectiveScheme)
```
The minimal concrete type realizing this interface is
```julia
ProjectiveScheme{CoeffRingType, RingType} <: AbsProjectiveScheme{CoeffRingType, RingType}
```


## Constructors

Besides `proj(S)` for some graded polynomial ring or a graded affine algebra `S`, we
provide the following constructors:
```julia
proj(S::MPolyDecRing)
proj(S::MPolyDecRing, I::MPolyIdeal{T}) where {T<:MPolyDecRingElem}
proj(I::MPolyIdeal{<:MPolyDecRingElem})
proj(Q::MPolyQuoRing{<:MPolyDecRingElem})
```
Subschemes defined by homogeneous ideals, ring elements, or lists of elements can be created
via the respective methods of the `subscheme(P::AbsProjectiveScheme, ...)` function.
Special constructors are provided for projective space itself via the function
`projective_space` and its various methods.
```@docs
projective_space(A::Ring, var_symb::Vector{VarName})
projective_space(A::Ring, r::Int; var_name::VarName=:s)
```

## Attributes
Besides those attributes already covered by the above general interface we have the following
(self-explanatory) ones for projective schemes over a field.
```julia
dim(P::AbsProjectiveScheme{<:Field})
hilbert_polynomial(P::AbsProjectiveScheme{<:Field})
degree(P::AbsProjectiveScheme{<:Field})
```

```@docs
arithmetic_genus(P::AbsProjectiveScheme{<:Field})
```

## Methods

To facilitate the interplay between an `AbsProjectiveScheme` and the affine charts of its
`covered_scheme` we provide the following methods:
```@docs
dehomogenization_map(X::AbsProjectiveScheme, U::AbsAffineScheme)
homogenization_map(P::AbsProjectiveScheme, U::AbsAffineScheme)
```

## Properties

Further properties of projective schemes:
```@docs
is_smooth(P::AbsProjectiveScheme{<:Any, <:MPolyQuoRing}; algorithm=:default)
```
```julia
is_empty(P::AbsProjectiveScheme{<:Field})
is_irreducible(P::AbsProjectiveScheme)
is_reduced(P::AbsProjectiveScheme)
is_geometrically_reduced(P::AbsProjectiveScheme{<:Field})
is_geometrically_irreducible(P::AbsProjectiveScheme{<:Field})
is_integral(X::AbsProjectiveScheme{<:Field})
is_geometrically_integral(X::AbsProjectiveScheme{<:Field})
```
