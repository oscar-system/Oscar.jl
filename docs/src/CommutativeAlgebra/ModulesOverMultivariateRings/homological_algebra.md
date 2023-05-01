```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

# Homological Algebra

Some OSCAR functions which are fundamental to homological algebra such as the `kernel` function
and basic functions for handling chain and cochain complexes have already been discussed
in previous sections. Building on these functions, we now introduce further OSCAR functionality
supporting computations in homological algebra.

## Presentations

```@docs
presentation(M::ModuleFP)
```

## Syzygies and Free Resolutions

```@docs
free_resolution(M::SubquoModule{<:MPolyRingElem}; 
    ordering::ModuleOrdering = default_ordering(M),
    length::Int=0, algorithm::Symbol=:fres
  )
```

## Homology

```@docs
homology(C::ComplexOfMorphisms{<:ModuleFP})
```

```@docs
homology(C::ComplexOfMorphisms{<:ModuleFP}, i::Int)
```

## Hom and Ext

```@docs
hom(M::ModuleFP, N::ModuleFP, alg::Symbol=:maps)
```

```@docs
element_to_homomorphism(f::ModuleFPElem)
```

```@docs
homomorphism_to_element(H::ModuleFP, phi::ModuleFPHom)
```

```@docs
ext(M::ModuleFP, N::ModuleFP, i::Int)
```

## Tensorproduct and Tor

```@docs
tensor_product(G::ModuleFP...; task::Symbol = :none)
```

```@docs
tor(M::ModuleFP, N::ModuleFP, i::Int)
```

## Fitting Ideals

```@docs
fitting_ideal(M::ModuleFP{T}, i::Int) where T <: MPolyRingElem
```

## Flatness

Checking flatness in OSCAR relies on characterizing flatness in terms of Fitting ideals.

```@docs
is_flat(M::ModuleFP{T}) where T <: MPolyRingElem
```

```@docs
non_flat_locus(M::ModuleFP{T}) where T <: MPolyRingElem
```

## Regular Sequence Test

```@docs
is_regular_sequence(V::Vector{T}, M::ModuleFP{T}) where T <: MPolyRingElem
```

## Koszul Complex

```@docs
koszul_matrix(V::Vector{T}, i::Int) where T <: MPolyRingElem
```

```@docs
koszul_complex(V::Vector{T}) where T <: MPolyRingElem
```

## Koszul Homology

```@docs
koszul_homology(V::Vector{T}, M::ModuleFP{T}, p::Int) where T <: MPolyRingElem
```

## Depth

The computation of depth in OSCAR relies on expressing depth in terms of  Koszul cohomology. 

```@docs
depth(I::MPolyIdeal{T}, M::ModuleFP{T}) where T <: MPolyRingElem
```






