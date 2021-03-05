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

## Ideal constructors

## Simple ideal Operations

### Powers of Ideal

```@docs
:^(I::MPolyIdeal, m::Int)
```
### Sum of Ideals

```@docs
:+(I::MPolyIdeal, J::MPolyIdeal)
```

### Product of Ideals

```@docs
:*(I::MPolyIdeal, J::MPolyIdeal)
```
###### Example

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
I = ideal(R, [x, y])
J = ideal(R, [z^2])
I^3
I+J
I*J
```

## Module Constructors

## Simple Module Operations

### Sum of Submodules

    N_1+N_2

### Products of Ideal and Submodule

    I*N

## Gröbner Bases

### Monomial Orders

### Normal Forms

### Computing Gröbner Bases

#### Gröbner Bases over fields

    groebner_basis(I)

    reduced_groebner_basis(I)

    fglm

    Gröbner walks

    Hilbert-driven

#### Leading Ideals

```@docs
leading_ideal(g::Array{T, 1}, args...) where { T <: MPolyElem }
```

###### Example

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
L = leading_ideal([x*y^2-3*x, x^3-14*y^5])
```  

```@docs
leading_ideal(I::MPolyIdeal)
```

###### Example

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
I = ideal(R,[x*y^2-3*x, x^3-14*y^5])
L = leading_ideal(I)
L = leading_ideal(I, :lex)
```  

#### Gröbner Bases over the integers

#### ....


## Ideals in Multivariate Polynomial Rings

### Introduction

Distinction: ideals <--> homogeneous ideals. To be discussed.

### Ideal Operations for Multivariate Polynomial Rings

#### Intersection of Ideals

```@docs
intersect(I::MPolyIdeal, J::MPolyIdeal)
```

###### Example

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
I = intersect(ideal(R, [x, y])^2, ideal(R, [y^2-x^3+x]))
```  

#### Ideal Quotients

```@docs
quotient(I::MPolyIdeal, J::MPolyIdeal)
```

###### Example

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
I = ideal(R, [x^4+x^2*y*z+y^3*z, y^4+x^3*z+x*y^2*z, x^3*y+x*y^3])
J = ideal(R,[x,y,z])^2
L = quotient(I,J)
I:J
```
#### Saturation

```@docs
saturation(I::MPolyIdeal, J::MPolyIdeal)
```

###### Example

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
I = ideal(R, [x^4+x^2*y*z+y^3*z, y^4+x^3*z+x*y^2*z, x^3*y+x*y^3])
J = ideal(R, [x,y,z])^2
L = saturation(I,J)
```


#### Elimination

```@docs
eliminate(I::MPolyIdeal, polys::Array{MPolyElem, 1})
```

```@docs
eliminate(I::MPolyIdeal, polys::AbstractArray{Int, 1})
```

###### Example

```@repl oscar
R, (t, x, y, z) = PolynomialRing(QQ, ["t", "x", "y", "z"])
I = ideal(R, [t-x, t^2-y, t-z^3])
A = [t]
AA = [1]
TC = eliminate(I,A)
TC = eliminate(I,AA)
```	
#### Homogenization and Dehomogenization

    homogenize(I,t)   CAVEAT: Als Ideal! Auch poly, vector, etc. Siehe M2.

    dehomogenize(I,t)

#### Decomposition

#### Radical

```@docs
radical(I::MPolyIdeal)
```

###### Example

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
I = intersect(ideal(R, [x, y])^2, ideal(R, [y^2-x^3+x]))
I = intersect(I, ideal(R, [x-y-1])^2)
RI = radical(I)

R, (a, b, c, d) = PolynomialRing(ZZ, ["a", "b", "c", "d"])
I = intersect(ideal(R, [9,a,b]), ideal(R, [3,c]))
I = intersect(I, ideal(R, [11,2a,7b]))
I = intersect(I, ideal(R, [13a^2,17b^4]))
I = intersect(I, ideal(R, [9c^5,6d^5]))
I = intersect(I, ideal(R, [17,a^15,b^15,c^15,d^15]))
RI = radical(I)
```

##### Primary Decomposition

```@docs
primary_decomposition(I::MPolyIdeal)
```

###### Example

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
I = intersect(ideal(R, [x, y])^2, ideal(R, [y^2-x^3+x]))
I = intersect(I, ideal(R, [x-y-1])^2)
L = primary_decomposition(I)

L = primary_decomposition(I, alg=:SY)

R, (a, b, c, d) = PolynomialRing(ZZ, ["a", "b", "c", "d"])
I = ideal(R, [1326*a^2*d^5, 1989*a^2*c^5, 102*b^4*d^5, 153*b^4*c^5,
663*a^2*c^5*d^5, 51*b^4*c^5*d^5, 78*a^2*d^15, 117*a^2*c^15,
78*a^15*d^5, 117*a^15*c^5, 6*a^2*b^4*d^15, 9*a^2*b^4*c^15,
39*a^2*c^5*d^15, 39*a^2*c^15*d^5, 6*a^2*b^15*d^5, 9*a^2*b^15*c^5,
6*a^15*b^4*d^5, 9*a^15*b^4*c^5, 39*a^15*c^5*d^5, 3*a^2*b^4*c^5*d^15,
3*a^2*b^4*c^15*d^5, 3*a^2*b^15*c^5*d^5, 3*a^15*b^4*c^5*d^5])
L = primary_decomposition(I)
```	

##### Minimal Associated Primes

```@docs
minimal_primes(I::MPolyIdeal)
```

###### Example

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
I = intersect(ideal(R, [x, y])^2, ideal(R, [y^2-x^3+x]))
I = intersect(I, ideal(R, [x-y-1])^2)
L = minimal_primes(I)

L = minimal_primes(I, alg=:charSets)

R, (a, b, c, d) = PolynomialRing(ZZ, ["a", "b", "c", "d"])
I = ideal(R, [1326*a^2*d^5, 1989*a^2*c^5, 102*b^4*d^5, 153*b^4*c^5,
663*a^2*c^5*d^5, 51*b^4*c^5*d^5, 78*a^2*d^15, 117*a^2*c^15,
78*a^15*d^5, 117*a^15*c^5, 6*a^2*b^4*d^15, 9*a^2*b^4*c^15,
39*a^2*c^5*d^15, 39*a^2*c^15*d^5, 6*a^2*b^15*d^5, 9*a^2*b^15*c^5,
6*a^15*b^4*d^5, 9*a^15*b^4*c^5, 39*a^15*c^5*d^5, 3*a^2*b^4*c^5*d^15,
3*a^2*b^4*c^15*d^5, 3*a^2*b^15*c^5*d^5, 3*a^15*b^4*c^5*d^5])
L = minimal_primes(I)
```
##### Weak Equidimensional Decomposition

```@docs
equidimensional_decomposition_weak(I::MPolyIdeal)
```
###### Example

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
I = intersect(ideal(R, [x, y])^2, ideal(R, [y^2-x^3+x]))
I = intersect(I, ideal(R, [x-y-1])^2)
L = equidimensional_decomposition_weak(I)
```

##### Equidimensional Decomposition of radical

```@docs
equidimensional_decomposition_radical(I::MPolyIdeal)
```

###### Example

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
I = intersect(ideal(R, [x, y])^2, ideal(R, [y^2-x^3+x]))
I = intersect(I, ideal(R, [x-y-1])^2)
L = equidimensional_decomposition_radical(I)
```

##### Equidimensional Hull

```@docs
equidimensional_hull(I::MPolyIdeal)
```

###### Example

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
I = intersect(ideal(R, [x, y])^2, ideal(R, [y^2-x^3+x]))
I = intersect(I, ideal(R, [x-y-1])^2)
L = equidimensional_hull(I)
  
R, (a, b, c, d) = PolynomialRing(ZZ, ["a", "b", "c", "d"])
I = ideal(R, [1326*a^2*d^5, 1989*a^2*c^5, 102*b^4*d^5, 153*b^4*c^5,
663*a^2*c^5*d^5, 51*b^4*c^5*d^5, 78*a^2*d^15, 117*a^2*c^15,
78*a^15*d^5, 117*a^15*c^5, 6*a^2*b^4*d^15, 9*a^2*b^4*c^15,
39*a^2*c^5*d^15, 39*a^2*c^15*d^5, 6*a^2*b^15*d^5, 9*a^2*b^15*c^5,
6*a^15*b^4*d^5, 9*a^15*b^4*c^5, 39*a^15*c^5*d^5, 3*a^2*b^4*c^5*d^15,
3*a^2*b^4*c^15*d^5, 3*a^2*b^15*c^5*d^5, 3*a^15*b^4*c^5*d^5])
L = equidimensional_hull(I)
```
##### Radical of the Equidimensional Hull

```@docs
equidimensional_hull_radical(I::MPolyIdeal)
```
###### Example

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
I = intersect(ideal(R, [x, y])^2, ideal(R, [y^2-x^3+x]))
I = intersect(I, ideal(R, [x-y-1])^2)
L = equidimensional_hull_radical(I)
```
##### Absolute Primary Decomposition

    absolute_primary_decomposition(I)               --->  absPrimdecGTZ   


### Tests on Ideals
  
    homogeneity test:        ishomogeneous(I)
    equality test:                 I==J
	ideal membership test: iscontained(f,I) oder    f \in I
    containment test:          issubset(I,J)    oder    I \subseteq J
    .....
    
### Data Associated to Ideals

    generators:            gens(I)                                (type list)
    dimension:             dim(I)                                 (type integer)
    codimension:           codim(I)                               (type integer)
    
### Data Associated to Homogeneous Ideals
    
                             min_base(I)   (or so)
     
    degree:                  degree(I)                             (type integer)
    Hilbert function:        hilbert_function(d,I)                 (type integer)
    Hilbert series:          hilbert_series(I)                     (type univariate rational function)  numerator, denominator
    reduced Hilbert series:  reduced_Hilbert_series(I)             (type univariate rational function)
    Hilbert polynomial:      hilbert_polynomial(I)                 (type univariate polynomial, direkt in Julia von hilbert_series(I))


## Modules over Multivariate Polynomial Rings

### Flatness

### Depth and Codimension

## Normalization of Rings

```@docs
normalize(Q::MPolyQuo)
```

###### Example

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
Q, _ = quo(R, ideal(R, [z - x^4, z - y^6]))
L = normalize(Q)
```

## Properties of Rings

 ### Normality

    isnormal(R)

### Cohen-Macaulayness

    iscohenmacaulay(R)

## Invariant Theory

### Invariant Theory of Finite Groups

### ...

