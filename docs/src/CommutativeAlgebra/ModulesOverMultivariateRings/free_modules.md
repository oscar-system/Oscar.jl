```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# [Free Modules](@id free_modules)

In this section, the expression *free module*  refers to a free module of finite rank
over a ring of type `MPolyRing`, `MPolyQuoRing`, `MPolyLocRing`, or `MPolyQuoLocRing`.
More concretely, given a ring $R$ of one of these types, the free $R$-modules considered are of
type $R^p$, where we think of $R^p$ as a free module with a given basis, namely the basis of
standard unit vectors. Accordingly, elements of free modules are represented by coordinate vectors,
and homomorphisms between free modules by matrices.

!!! note
    By convention, vectors are row vectors, and matrices operate by multiplication on the right.

## Types

All OSCAR types for the modules considered here belong to the
abstract type `ModuleFP{T}`, where `T` is the element type of the underlying ring.
Graded or not, the free modules belong to the abstract subtype `AbstractFreeMod{T} <: ModuleFP{T}`,
they are modeled as objects of the concrete type `FreeMod{T} <: AbstractFreeMod{T}`.

!!! note
    Canonical maps such us the canonical projection onto a quotient module arise in many 
    constructions in commutative algebra. The `FreeMod` type is designed so that it allows
    for the caching of such maps when executing functions. The `direct_sum` function discussed
    in this section provides an example.

## Constructors

```@docs
free_module(R::MPolyRing, n::Int, name::VarName = :e; cached::Bool = false)
```

Over graded multivariate polynomial rings and their quotients,  there are two basic ways of
creating graded free modules: While the `grade` function allows one to create a graded free module
by assigning a grading to a free module already constructed, the `graded_free_module` function is
meant to create a graded free module all at once.

```@docs
grade(F::FreeMod, W::Vector{FinGenAbGroupElem})
```

```@docs
grade(F::FreeMod, W::Vector{<:Vector{<:IntegerUnion}})
```

```@docs
graded_free_module(R::Ring, p::Int, W::Vector{FinGenAbGroupElem}=[grading_group(R)[0] for i in 1:p], name::String="e")
```

```@docs
graded_free_module(R::Ring, W::Vector{<:Vector{<:IntegerUnion}}, name::String="e")
```

## Data Associated to Free Modules

If `F` is a free `R`-module, then

- `base_ring(F)` refers to `R`,
- `basis(F)`, `gens(F)` to the basis vectors of `F`, 
- `rank(F)`, `number_of_generators(F)` / `ngens(F)`, `dim(F)` to the number of these vectors, and
- `F[i]`, `basis(F, i)`, `gen(F, i)` to the `i`-th such vector.

###### Examples

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> F = free_module(R, 3);

julia> basis(F)
3-element Vector{FreeModElem{QQMPolyRingElem}}:
 e[1]
 e[2]
 e[3]

julia> rank(F)
3
```

In the graded case, we also have:

```@docs
 grading_group(F::FreeMod)
```

```@docs
degrees_of_generators(F::FreeMod)
```

## Elements of Free Modules

All OSCAR types for elements of the modules considered here belong
to the abstract type `ModuleElemFP{T}`, where `T` is the element type of the underlying ring.
The free modules belong to the abstract subtype `AbstractFreeModElem{T} <: ModuleFPElem{T}`.
They are modeled as objects of the concrete type `FreeModElem{T} <: AbstractFreeModElem{T}`
which implements an element $f$ of a free module $F$ as a sparse row, that is, as an object of
type `SRow{T}`. This object specifies the coordinates of $f$ with respect to the basis of standard
unit vectors of $F$. To create an element, enter its coordinates as a sparse row or a vector: 

```@julia
(F::FreeMod{T})(c::SRow{T}) where T
```

```@julia
(F::FreeMod{T})(c::Vector{T}) where T
```

Alternatively, directly write the element as a linear combination of basis vectors of $F$:
 
##### Examples

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> F = free_module(R, 3);

julia> f = F(sparse_row(R, [(1,x),(3,y)]))
x*e[1] + y*e[3]

julia> g = F([x, zero(R), y])
x*e[1] + y*e[3]

julia> h = x*F[1] + y*F[3]
x*e[1] + y*e[3]

julia> f == g == h
true
```

Given an element `f`  of a free module `F` over a multivariate polynomial ring with element type `T`,
- `parent(f)` refers to `F`, and
- `coordinates(f)` to the coordinate vector of `f`, returned as an object of type `SRow{T}`.

##### Examples

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> F = free_module(R, 3);

julia> f = x*F[1] + y*F[3]
x*e[1] + y*e[3]

julia> parent(f)
Free module of rank 3 over R

julia> coordinates(f)
Sparse row with positions [1, 3] and values QQMPolyRingElem[x, y]

```

The zero element of a free module is obtained as follows:

```@docs
zero(F::AbstractFreeMod)
```

Whether a given element of a free module is zero can be tested as follows:

```@docs
is_zero(f::AbstractFreeModElem)
```

In the graded case, we additionally have:

```@docs
is_homogeneous(f::FreeModElem)
```

```@docs
degree(f::FreeModElem{T}) where {T<:Union{<:MPolyDecRingElem, <:MPolyQuoRingElem{<:MPolyDecRingElem}}}
```


## Tests on Free Modules

The tests [`is_graded`](@ref), [`is_standard_graded`](@ref), [`is_z_graded`](@ref),
and [`is_zm_graded`](@ref) carry over analogously to free modules. They return `true` if the
corresponding property is satisfied, and `false` otherwise. In addition, we have:

```@docs
==(F::FreeMod, G::FreeMod)
```

```@docs
is_isomorphic(F::FreeMod, G::FreeMod)
```

```@docs
is_zero(F::AbstractFreeMod)
```

## Homomorphisms from Free Modules

All OSCAR types for homomorphisms of the modules considered here belong
to the abstract type `ModuleFPHom{T1, T2}`, where `T1` and `T2` are the types of domain and codomain respectively.
A homomorphism $F\to M$ from a free module $F$ is determined by specifying the images
of the basis vectors of $F$ in $M$. For such homomorphisms, OSCAR provides the concrete type
`FreeModuleHom{T1, T2} <: ModuleFPHom{T1, T2}` as well as the following constructors:

```@docs
hom(F::FreeMod, M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}) where T 
```

```@docs
hom(F::FreeMod, M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}, h::RingMapType) where {T, RingMapType}
```

Given a homomorphism of type `FreeModuleHom`, a matrix representing it
is recovered by the following function:

```@docs
matrix(a::FreeModuleHom)
```

The domain and codomain of a homomorphism `a`  of type `FreeModuleHom` can be
recovered by entering `domain(a)` and `codomain(a)`, respectively.

The functions below test whether a homomorphism of type
`FreeModuleHom` is graded and homogeneous, respectively.

```@docs
is_graded(a::FreeModuleHom)
```

```@docs
is_homogeneous(a::FreeModuleHom)
```

In the graded case, we additionally have:

```@docs
degree(a::FreeModuleHom)
```

```@docs
grading_group(a::FreeModuleHom)
```


