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

$a:R^m ⟶R^p \;\text{ and }\; b:R^n ⟶R^p$

are two homomorphisms of free $R$-modules with the same codomain. We then refer to
- the module $M$ as the *subquotient defined by $a$ and $b$*,
- the codomain $R^p$ as the *ambient free module* of $M$,
- the images of the canonical basis vectors of $R^m$ in $M$ as the *generators* of $M$, and
- the images of the canonical basis vectors of $R^n$ in $R^p$ as the *relations* of $M$.

Alternatively, we speak of the *subquotient of* $\text{im } a$ *by* $\text{im } b$ or the
*subquotient defined by $A$ and $B$*, where $A$ and $B$ are the matrices representing
$a$ and $b$, respectively. 

!!! note
    Recall from the section on free modules that by a free $R$-module we mean a free
    module of type $R^n$ , where we think of $R^n$ as a free module with a given
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
Its concrete descendant `SubQuoElem{T}` implements an element $f$ of a subquotient
$M$ over the ring $R$ as a sparse row, that is, as an object of type `SRow{T}`, which
specifies the coefficients of an $R$-linear combination of the generators of $M$
defining $f$. To create an element, simply write it as an $R$-linear combination
of generators of $M$, or enter its coefficients as an object of type
`SRow{T}` or `Vector{T}`: 

```@julia
(M::SubQuo{T})(c::SRow{T}) where T
```

```@julia
(M::SubQuo{T})(c::Vector{T}) where T
```

##### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
A = R[x; y]
B = R[x^2; x*y; y^2; z^4]
M = SubQuo(A, B)
f = z*M[1] + M[2]
g = M(sparse_row(R, [(1,z),(2,one(R))]))
h = M([z, one(R)])
f == g == h
```



## Gröbner Bases



## Tests on Subqotients

### Containment of Subqotients


### Equality of Subqotients


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







