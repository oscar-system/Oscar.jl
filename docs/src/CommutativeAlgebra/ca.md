```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["integer.md"]
```

# Commutative Algebra

## Introduction

This part of OSCAR is under development with regard to providing both the functionality and the documentation. At current state, the documentation serves as a guide on what to do and will be constantly updated.

General textbooks offering details on theory and algorithms include: 
- [GP07](@cite)
- [DL06](@cite)
- [DP13](@cite)

## Ideal Constructors for Multivariate Polynomial Rings

```@docs
ideal(g::Array{T, 1}) where {T <: MPolyElem}
```
```@docs
ideal(Rx::MPolyRing, g::Array{<:Any, 1})
```

## Gröbner Bases

### Monomial Orders

### Normal Forms

### Computing Gröbner Bases

```@docs
groebner_basis(I::MPolyIdeal)
```

```@docs
groebner_basis(I::MPolyIdeal, ord::Symbol; complete_reduction::Bool=false)
```

#### Gröbner Bases with transformation matrix

```@docs
groebner_basis_with_transformation_matrix(I::MPolyIdeal; ord::Symbol = :degrevlex, complete_reduction::Bool=false)
```

    fglm

    Gröbner walks

    Hilbert-driven

#### Leading Ideals

```@docs
leading_ideal(g::Array{T, 1}, args...) where { T <: MPolyElem }
```

```@docs
leading_ideal(I::MPolyIdeal)
```

#### Gröbner Bases over the integers

#### ....


### Syzygies

#### Generators of syzygies

```@docs
syzygy_generators(a::Array{<:MPolyElem, 1})
```

## Ideals in Multivariate Polynomial Rings

### Introduction

Distinction: ideals <--> homogeneous ideals. To be discussed.

### Ideal Operations for Multivariate Polynomial Rings

#### Simple ideal Operations

##### Powers of Ideal

```@docs
:^(I::MPolyIdeal, m::Int)
```
##### Sum of Ideals

```@docs
:+(I::MPolyIdeal, J::MPolyIdeal)
```

##### Product of Ideals

```@docs
:*(I::MPolyIdeal, J::MPolyIdeal)
```

#### Intersection of Ideals

```@docs
intersect(I::MPolyIdeal, J::MPolyIdeal)
```

#### Ideal Quotients

```@docs
quotient(I::MPolyIdeal, J::MPolyIdeal)
```

#### Saturation

```@docs
saturation(I::MPolyIdeal, J::MPolyIdeal)
```

#### Elimination

```@docs
eliminate(I::MPolyIdeal, polys::Array{MPolyElem, 1})
```

```@docs
eliminate(I::MPolyIdeal, polys::AbstractArray{Int, 1})
```

#### Homogenization and Dehomogenization

    homogenize(I,t)   CAVEAT: Als Ideal! Auch poly, vector, etc. Siehe M2.

    dehomogenize(I,t)

#### Decomposition

#### Radical

```@docs
radical(I::MPolyIdeal)
```

##### Primary Decomposition

```@docs
primary_decomposition(I::MPolyIdeal)
```

##### Minimal Associated Primes

```@docs
minimal_primes(I::MPolyIdeal)
```

##### Weak Equidimensional Decomposition

```@docs
equidimensional_decomposition_weak(I::MPolyIdeal)
```

##### Equidimensional Decomposition of radical

```@docs
equidimensional_decomposition_radical(I::MPolyIdeal)
```

##### Equidimensional Hull

```@docs
equidimensional_hull(I::MPolyIdeal)
```

##### Radical of the Equidimensional Hull

```@docs
equidimensional_hull_radical(I::MPolyIdeal)
```

##### Absolute Primary Decomposition

    absolute_primary_decomposition(I)               --->  absPrimdecGTZ   


### Tests on Ideals

#### Equality of Ideals

```@docs
:(==)(I::MPolyIdeal, J::MPolyIdeal)
```

#### Containment of Ideals

```@docs
issubset(I::MPolyIdeal, J::MPolyIdeal)
```

#### Ideal Membership

```@docs
ideal_membership(f::MPolyElem, I::MPolyIdeal)
```

#### Radical Membership

```@docs
radical_membership(f::MPolyElem, I::MPolyIdeal)
```

#### Primality Test

```@docs
isprime(I::MPolyIdeal)
```

#### Primary Test

```@docs
isprimary(I::MPolyIdeal)
```
  
    homogeneity test:        ishomogeneous(I)
    
    .....
    
### Data Associated to Ideals

```@docs
base_ring(I::MPolyIdeal)
```

#### Number of Generators

```@docs
ngens(I::MPolyIdeal)
```

#### Generators

```@docs
gens(I::MPolyIdeal)
```

#### Dimension

```@docs
dim(I::MPolyIdeal)
```

#### Codimension

```@docs
codim(I::MPolyIdeal)
```
    
### Data Associated to Homogeneous Ideals
    
                             min_base(I)   (or so)
     
    degree:                  degree(I)                             (type integer)
    Hilbert function:        hilbert_function(d,I)                 (type integer)
    Hilbert series:          hilbert_series(I)                     (type univariate rational function)  numerator, denominator
    reduced Hilbert series:  reduced_Hilbert_series(I)             (type univariate rational function)
    Hilbert polynomial:      hilbert_polynomial(I)                 (type univariate polynomial, direkt in Julia von hilbert_series(I))


## Modules over Multivariate Polynomial Rings

### Module Constructors

### Simple Module Operations

#### Sum of Submodules

    N_1+N_2

#### Products of Ideal and Submodule

    I*N

### Flatness

### Depth and Codimension

## Affine Algebras

An affine algebra over a field $K$, or an affine $K$-algebra, is the quotient $A=R/I$ of a multivariate
polynomial ring $R$ over $K$ modulo a radical ideal $I$.

### Constructor

```@docs
quo(R::MPolyRing, I::MPolyIdeal)
```
###### Example

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
A, _ = quo(R, ideal(R, [x^2-y^3, x-y]))
A, f = quo(R, ideal(R, [x^2-y^3, x-y]))
f
```

### Data Associated to Affine Algebras

#### Basic Data

If `B = S/J` is an affine algebra, then `B.R` refers to `S` and `A.I` to `J`.

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
B, _ = quo(R, ideal(R, [y-x^2, x-z^3]))
B.R
B.I
```

Furthermore, `gens(B)` refers to the generators of `B` and  `ngens(B)` to their number.

```@repl oscar
gens(B)
ngens(B)
```

#### Dimension

```@docs
dim(A::MPolyQuo)
```

###### Example

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
A, _ = quo(R, ideal(R, [y-x^2, x-z^3]))
dim(A)
```

### Ideals in Affine Algebras


### Tests on Affine Algebras

#### Subalgebra Membership

```@docs
subalgebra_membership(A::MPolyQuo, f::MPolyQuoElem, v::Vector{T}) where T <: MPolyQuoElem
```

#### Reducedness Test

```@docs
isreduced(A::MPolyQuo)
```

###### Example

```@repl oscar
R, (x,) = PolynomialRing(QQ, ["x"])
A, _ = quo(R, ideal(R, [x^4]))
isreduced(A)
```

#### Normality Test

```@docs
isnormal(A::MPolyQuo)
```

#### Cohen-Macaulayness Test

iscohenmacaulay(R)

###### Example

### Homomorphisms of Affine Algebras

#### Constructors

```@docs
hom(D::MPolyQuo, C::MPolyQuo, V::Vector)
```
###### Example

#### Data Associated to Homomorphisms of Affine Algebras

```@docs
domain(F::MPolyQuoHom)
```

```@docs
codomain(F::MPolyQuoHom)
```
#### Tests on Homomorphisms of Affine Algebras

```@docs
isinjective(F::MPolyQuoHom)
```

```@docs
issurjective(F::MPolyQuoHom)
```

```@docs
isbijective(F::MPolyQuoHom)
```

```@docs
isfinite(F::MPolyQuoHom)
```

### Noether Normalization

```@docs
noether_normalization(A::MPolyQuo)
```

###### Example

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
A, _ = quo(R, ideal(R, [x*y, x*z]))
L = noether_normalization(A);
L[1]
L[2]
L[3]
```

### Normalization of Rings

```@docs
normalize(A::MPolyQuo)
```

###### Example

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
A, _ = quo(R, ideal(R, [(x^2-y^3)*(x^2+y^2)*x]))
L = normalize(A)

R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
A, _ = quo(R, ideal(R, [z^3-x*y^4]))
L = normalize(A)
```



#### 

## Invariant Theory

### Invariant Theory of Finite Groups

### ...
