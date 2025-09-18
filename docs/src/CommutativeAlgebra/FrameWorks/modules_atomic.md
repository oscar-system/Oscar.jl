```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Atomic Functions for Subquotient Module Framework

This section documents the atomic functions on the level of submodules of free modules which are the minimal
requirement for the subquotient module framework and the homological methods built on top of it to function properly. Developers should implement these functions specifically for modules over their ring type under consideration. We assume in the following that your custom ring type is `MyRing` and
the element type is `MyRingElem`. 

---

## Membership Test

```@docs
in_atomic(a::FreeModElem, M::SubModuleOfFreeModule)
```

**Implementation Guidance:**

The signature implemented should look like this:

```julia
function in_atomic(a::FreeModElem{T}, M::SubModuleOfFreeModule) where {T<:MyRingElem}
```

---

## Coordinates Computation

```@docs
coordinates_atomic(a::FreeModElem, M::SubModuleOfFreeModule; task::Symbol = :auto)
```

**Implementation Guidance:**

The signature implemented should look like this:

```julia
function coordinates_atomic(a::FreeModElem{T}, M::SubModuleOfFreeModule; task::Symbol=:auto) where {T<:MyRingElem}
```

---

## Kernel Computation (`kernel_atomic`)

```@docs
kernel_atomic(h::FreeModuleHom{<:FreeMod, <:FreeMod})
```

**Implementation Guidance:**

The signature implemented should look like this:

```julia
function kernel_atomic(h::FreeModuleHom{<:FreeMod{T}, <:FreeMod{T}, Nothing}) where {T<:MyRingElem}
```

---

## Pruning a Module, also supplying the Isomorphism (`prune_with_map_atomic`)

```@docs
prune_with_map_atomic(M::ModuleFP)
```

**Implementation Guidance:**

There is a generic fallback for this function (which is not just the identity), but since this function is used to simplify presentations
and hence is crucial for performance one should always try to make it as efficient as possible 
(both with respect to its performance and size of its answer).

The signature implemented should look like this:

```julia
function prune_with_map_atomic(M::SubquoModule{T}) where {T<:MyRingElem}
```

---

## Optional Functions

Various functions can optionally be implemented for performance benefits or to handle special cases. Examples:

- `free_resolutions`: Provide a specific implementation of free resolutions (potentially faster than iterated kernel computation).
- `size`: Compute the cardinality in case the module is finite.
- `is_finite`: Check finiteness as a set.
- An iterator over elements or basis elements, if applicable.

For example, the signature implemented for free resolutions should look like this:

```julia
function free_resolutions(M::SubquoModule{T}; length::Int = 0, algorithm::Symbol = :auto) where {T<:MyRingElem}
```

## Alternative: Implement a Gröbner basis framework

Where applicable, one can also opt to implement a Gröbner basis framework situated one abstraction layer below the submodule framework. This is at the current stage more intricate, please do not hesitate to get into contact.

