```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["galois.md"]
```

# Galois Theory

Let `K` be a finite (separable) extension of `k`. Then, in contrast to most of 
the literature we distinguish two concepts

 - the *automorphism group*
 - the *Galois group*
 
The automorphism group deals with the actual automorphism of `K` fixing `k`
and thus is, in general trivial. Access to via two constructions:

 - a list of all automorphisms (usually on the identity)
 - the group of automorphisms, returned as an abstact group and
   a map linking group elements to actual automorphisms

On the other hand, the Galois group is isomorphic to the automorphism
group of the normal closure and is explicitly given as a group of
permutations of the roots  of the defining polynomial. Thus even
in the case of `K` over `k` being normal, elements of the
Galois group do not immediately give automorphisms at all.

Currently, the computation of Galois groups is possible for
 
  - `K` a simple extension of the rationals (`AnticNumberField`)
  - `K` a simple extension of an `AnticNumberField`
  - `K` a finite extension of the rational function field over the
     rationals
  - `f` a polynomial over the rationals, an `AnticNumberField` or 
     the univarite function field over the rationals (here explicit
     permutations of the roots in a suitable splitting field are returned).
     Futhermore, the the case of `f` over the rational function field, the
     monodromy can be computed, ie. the automorphism group over the
     complex numbers.

Independently of the Galois group, subfields, that is intermediate fields
between `K` and `k` can be computed as well.

## Automorphism Group

The automorphisms are computed using various specialised factoring
algorithms: lifting the roots of the defining polynomial in the
given field modulo suitable prime ideal powers and
recovering the true roots from this information.

```@docs
automorphisms(::NumField)
automorphism_group(K::AnticNumberField)
automorphism_group(L::NumField, K::NumField)
absolute_automorphism_group(L::NumField)
```

```@docs
fixed_field(K::SimpleNumField, sigma::T; simplify::Bool = true) where {T <: Hecke.NumFieldMor}
```

## Subfields

```@docs
subfields(K::SimpleNumField; degree::Int = -1)
Hecke.principal_subfields(K::SimpleNumField)
subfields(FF::Generic.FunctionField{fmpq})
```

By setting `set_verbose_level(:Subfields, n::Int)` to 1 or 2
information about the progress can be obtained.

## Galois Group

The computation of Galois groups follows Stauduhars algorithm with many
improvements, see ... for an overview.

The entrire computation can also be thought of finding a describtion of the
splitting field of the polynomial. In fact, the information returned
can be used to verify any algebraic identity between the roots, and
find explicit subfields of the splitting field as well.

Information about the progress is available via
 
 - `:GaloisGroup`
 - `:Invariants`

```@docs
galois_group(K::AnticNumberField, extra::Int = 5; useSubfields::Bool = true, pStart::Int = 2*degree(K), prime::Int = 0)
galois_group(f::PolyElem{<:FieldElem})
```

The information returned consists always of a group `G` and a `GaloisCtx`: `C`.
Jointly, they can be used to further work with the information:

```@docs
roots(C::Oscar.GaloisGrp.GaloisCtx{Hecke.qAdicRootCtx}, pr::Int)
Oscar.GaloisGrp.upper_bound()
Oscar.GaloisGrp.isinteger()
Oscar.GaloisGrp.resolvent(C::Oscar.GaloisGrp.GaloisCtx, G::PermGroup, U::PermGroup, extra::Int = 5)
```

```@docs
galois_quotient(C::Oscar.GaloisGrp.GaloisCtx, Q::PermGroup)
galois_quotient(C::Oscar.GaloisGrp.GaloisCtx, d::Int)
galois_quotient(C::Oscar.GaloisGrp.GaloisCtx, d::Int, n::Int)
galois_quotient(f::PolyElem, p::Vector{Int})
fixed_field(GC::Oscar.GaloisGrp.GaloisCtx, U::PermGroup, extra::Int = 5)
minpoly(C::Oscar.GaloisGrp.GaloisCtx, I, extra::Int = 5)
```

```@docs
Oscar.GaloisGrp.cauchy_ideal(f::PolyElem{<:FieldElem})
Oscar.GaloisGrp.galois_ideal(C::Oscar.GaloisGrp.GaloisCtx, extra::Int = 5)
```
