```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["free_modules.md"]
```

# Free Modules

In what follows, the expression *free module*  refers to a free module of finite rank
over a multivariate polynomial ring. More concretely, given a multivariate polynomial ring $R$, 
the free $R$-modules considered are of type $R^p$, where we think of $R^p$ as a
free module with a given basis, namely the basis of standard unit vectors.
Accordingly, elements of free modules are represented by coordinate vectors,
and homomorphisms between free modules by matrices.

!!! note
    By convention, vectors are row vectors, and matrices operate by multiplication on the right.

## Types

All OSCAR types for finitely presented modules over multivariate polynomial rings belong to
the abstract type `ModuleFP{T}`, where `T` is the element type of the polynomial ring.
For free modules, OSCAR provides the abstract subtype `AbstractFreeMod{T} <: ModuleFP{T}`
and its concrete descendant `FreeMod{T <: RingElem}`.

!!! note
    Canonical maps such us the canonical projection onto a quotient module arise in many 
    constructions in commutative algebra. The `FreeMod` type is designed so that it allows
    for the caching of such maps when executing functions. The `direct_sum` function discussed
    in this section provides an example.

## Constructor

```@docs
free_module(R::MPolyRing, n::Int, name::String = "e"; cached::Bool = false)
```

## Data Associated to Free Modules

If `F` is a free `R`-module, then

- `base_ring(F)` refers to `R`,
- `basis(F)`, `gens(F)` to the basis vectors of `F`, 
- `rank(F)`, `ngens(F)`, `dim(F)` to the number of these vectors, and
- `F[i]`, `basis(F, i)`, `gen(F, i)` to the `i`-th such vector.

###### Examples

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
F = free_module(R, 3)
basis(F)
rank(F)
```

## Elements of Free Modules

All OSCAR types for elements of finitely presented modules over multivariate polynomial rings belong
to the abstract type `ModuleElemFP{T}`, where `T` is the element type of the polynomial ring.
For elements of free modules, there are the abstract subtype `AbstractFreeModElem{T} <: ModuleFPElem{T}` and its concrete
descendant `FreeModElem{T}` which implements an element $f$ of a free module $F$ as a sparse row,
that is, as an object of type `SRow{T}`. This object specifies the coordinates of $f$ with respect to
the basis of standard unit vectors of $F$. To create an element, enter its coordinates as a sparse row or a vector: 


```@julia
(F::FreeMod{T})(c::SRow{T}) where T
```

```@julia
(F::FreeMod{T})(c::Vector{T}) where T
```

Alternatively, directly write the element as a linear combination of basis vectors of $F$:
 
##### Examples

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
F = free_module(R, 3)
f = F(sparse_row(R, [(1,x),(3,y)]))
g = F( [x, zero(R), y])
h = x*F[1] + y*F[3]
f == g == h
```

Given an element `f`  of a free module `F` over a multivariate polynomial ring with element type `T`,
- `parent(f)` refers to `F`, and
- `coefficients(f)` to the coordinate vector of `f`, returned as an object of type `SRow{T}`.

##### Examples

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
F = free_module(R, 3)
f = x*F[1] + y*F[3]
parent(f)
coefficients(f)
```

The zero element of a free module is obtained as follows:

```@docs
zero(F::AbstractFreeMod)
```

Whether a given element of a free module is zero can be tested as follows:

```@docs
iszero(f::AbstractFreeModElem)
```

## Tests on Free Modules

```@docs
==(F::FreeMod, G::FreeMod)
```

```@docs
iszero(F::AbstractFreeMod)
```

## Homomorphisms from Free Modules

A homomorphism $F\rightarrow M$ from a free module $F$ is determined by specifying the images
of the basis vectors of $F$ in $M$. In OSCAR, such homomorphisms have type `FreeModuleHom{T1, T2}`, where
`T1` and `T2` are the element types of the domain and codomain, respectively. They are created
by using one of the following constructors:

```@docs
hom(F::FreeMod, M::ModuleFP, V::Vector{<:ModuleFPElem}) 
```

Given a homomorphism of type `FreeModuleHom`, a matrix `A` representing it
is recovered by the following function:

```@docs
matrix(a::FreeModuleHom)
```

##### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
F = free_module(R, 3)
G = free_module(R, 2)
V = [y*G[1], x*G[1]+y*G[2], z*G[2]]
a = hom(F, G, V)
A = matrix(a)
a(F[2])
```

The domain and codomain of a homomorphism `a`  of type `FreeModuleHom` can be
recovered by entering `domain(a)` and `codomain(a)`, respectively.



