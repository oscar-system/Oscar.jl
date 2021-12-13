```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["modules.md"]
```

# Modules Over Polynomial Rings

In this section, the term module will refer to a finitely presented module over a multivariate polynomial ring.
In OSCAR, the most general way of implementing such a module is that of a *subquotient*, that is,
as a submodule of a quotient of a free module. Explicitly, a *subquotient* $M$ over the ring $R$ is a module of type

$M = (\text{im } a + \text{im } b)/\text{im } b,$

where

$a:R^s ⟶R^p \;\text{ and }\; b:R^t ⟶R^p$

are two homomorphisms of free $R$-modules with the same codomain. We then refer to
- the module $M$ as the *subquotient defined by $a$ and $b$*,
- the codomain $R^p$ as the *ambient free module* of $M$,
- the images of the canonical basis vectors of $R^s$ in $M$ as the *generators* of $M$, and
- the images of the canonical basis vectors of $R^t$ in $R^p$ as the *relations* of $M$.

Alternatively, we speak of the *subquotient of* $\text{im } a$ *by* $\text{im } b$ or the
*subquotient defined by $A$ and $B$*, where $A$ and $B$ are the matrices representing
$a$ and $b$, respectively. 

!!! note
    Recall from the section on free modules that by a free $R$-module we mean a free
    module of type $R^p$ , where we think of $R^p$ as a free module with a given
	basis, namely the basis of standard unit vectors. Accordingly, elements of free modules
	are represented by coordinate vectors, and homomorphisms between free modules by
	matrices. Here, by convention, vectors are row vectors, and matrices operate by
	multiplication on the right.

## Types

For subquotients, OSCAR provides the abstract type `AbstractSubQuo{T}` and its concrete
descendant `SubQuo{T}`. The abstract supertype for all finitely presented modules over
commutative rings is `ModuleFP{T}`.

## Constructors

```@docs
subquotient(a::FreeModuleHom{T}, b::FreeModuleHom{T}) where T
```

## Data Associated to Subqotients

If `M` is a subquotient with ambient free `R`-module `F`, then

- `base_ring(M)` refers to `R`,
- `ambient_free_module(M)` to `F`,
- `gens(M)` to the generators of `M`, 
- `ngens(M)` to the number of these generators, 
- `gen(M, i)` to the `i`th such generator, and
- `rels(M)` to the relations of `M`.

##### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
A = R[x; y]
B = R[x^2; x*y; y^2; z^4]
M = SubQuo(A, B)
base_ring(M)
ambient_free_module(M)
gens(M)
ngens(M)
gen(M, 2)
rels(M)
```

## Elements of Subqotients

The abstract type for elements of subquotients in OSCAR is `AbstractSubQuoElem{T}`.
Its concrete descendant `SubQuoElem{T}` implements an element $m$ of a subquotient
$M$ over the ring $R$ as a sparse row, that is, as an object of type `SRow{T}`.
This object specifies the coefficients of an $R$-linear combination of the generators of $M$
which gives $m$. To create an element, enter the coefficients as an object of type `SRow{T}` or `Vector{T}`: 

```@julia
(M::SubQuo{T})(c::SRow{T}) where T
```

```@julia
(M::SubQuo{T})(c::Vector{T}) where T
```

Alternatively, directly write the element as an $R$-linear combination of generators of $M$.

##### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
A = R[x; y]
B = R[x^2; x*y; y^2; z^4]
M = SubQuo(A, B)
m = M(sparse_row(R, [(1,z),(2,one(R))]))
n = M([z, one(R)])
o = z*M[1] + M[2]
m == n == o
```

Given an element `m`  of a subquotient `M`,
- `parent(m)` refers to `M`, and
- `coefficients(m)` to  an object of type `SRow{T}` specifying the coefficients of an $R$-linear combination of the generators of $M$ which gives $m$. 

##### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);
A = R[x; y]
B = R[x^2; x*y; y^2; z^4]
M = SubQuo(A, B);
m = z*M[1] + M[2]
parent(m)
coefficients(m)
```

The zero element of a subquotient is obtained as follows:

```@docs
zero(M::SubQuo)
```

Whether a given element of a subquotient is zero can be tested as follows:

```@docs
iszero(m::SubQuoElem)
```

## Tests on Subqotients



## Basic Operations on Subqotients


## Submodules and Quotients



## Homomorphisms of Subqotients



## Operations on Homomorphisms of Subqotients



## Subquotients Related to Homomorphisms

### Kernel



### Image



### Cokernel


### Homology




## Presentations


## Syzygies and Free Resolutions



## Hom and Ext



## Tensorproduct and Tor







