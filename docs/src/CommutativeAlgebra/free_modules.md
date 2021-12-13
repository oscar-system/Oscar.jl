```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["free_modules.md"]
```

# Free Modules Over Commutative Rings

In this section, the expression *free module*  refers to a free module of finite rank
over a commutative ring. More concretely, given a commutative ring $R$, 
the free $R$-modules we consider are of type $R^p$, where we think of $R^p$ as a
free module with a given basis, namely the basis of standard unit vectors.
Accordingly, elements of free modules are represented by coordinate vectors,
and homomorphisms between free modules by matrices.

!!! note
    By convention, vectors are row vectors, and matrices operate by multiplication on the right.

## Types

For free modules, OSCAR provides the abstract type `AbstractFreeMod{T}` and its concrete descendant `FreeMod{T <: RingElem}`.
The abstract supertype for all finitely presented modules over commutative rings is `ModuleFP{T}`.


## Constructor

```@docs
free_module(R::Ring, n::Int, name::String = "e"; cached::Bool = false)
```

## Data Associated to Free Modules

If `F` is a free `R`-module, then

- `base_ring(F)` refers to `R`,
- `basis(F)`, `gens(F)` to the basis vectors of `F`, 
- `rank(F)`, `ngens(F)`, `dim(F)` to the number of these vectors, and
- `F[i]`, `basis(F, i)`, `gen(F, i)` to the `i`th such vector.

###### Examples

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
F = free_module(R, 3)
basis(F)
rank(F)
```

## Elements of Free Modules

The abstract type for elements of free modules in OSCAR is `AbstractFreeModElem{T}`. Its concrete
descendant `FreeModElem{T}` implements an element $f$ of a free module $F$ as a sparse row,
that is, as an object of type `SRow{T}`. This object specifies the coordinates of $f$ with respect to
the basis of standard unit vectors of $F$. To create an element, enter its coordinates as an object
of type `SRow{T}` or `Vector{T}`: 


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

Given an element `f`  of a free module `F`,
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

## Tests on Free Modules

```@docs
==(F::FreeMod, G::FreeMod)
```

```@docs
iszero(F::AbstractFreeMod)
```

## Homomorphisms From Free Modules

A homomorphism from a free module `F` is determined by specifying the images
of the basis vectors of `F`. In OSCAR, such homomorphisms have type `FreeModuleHom{T1, T2}`, where
`T1` and `T2` are the element types of the domain and codomain, respectively. They can be
constructed as follows:

```@docs
hom(F::FreeMod, G::ModuleFP, V::Vector)
```

If both `F` and `G`  are free modules, a  homomorphism `F` $\to$ `G` can be created
by entering the matrix representing it:

```@docs
hom(F::FreeMod{T}, G::ModuleFP{T}, A::MatElem{T}) where T
```

The matrix of a homomorphism `a` between free modules can be recovered by entering `matrix(a)`.

##### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
F = free_module(R, 3)
G = free_module(R, 2)
V = [y*G[1], x*G[1]+y*G[2], z*G[2]]
a = hom(F, G, V)
A = matrix(a)
```

If  `a`  is any homomorphism of type `FreeModuleHom`, its domain and codomain can be
recovered by entering `domain(a)` and `codomain(a)`, respectively.

## Operations on Free Modules

```@docs
direct_sum(F::FreeMod{T}...; task::Symbol = :none) where T
```

```@docs
direct_product(F::FreeMod{T}...; task::Symbol = :none) where T
```

##### Examples

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
F = free_module(R, 2)
G = free_module(R, 3, "f")
L = direct_product(F, G, task = :both)
basis(L[1])
L[3][2]
```

