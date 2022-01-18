```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["rings.md"]
```

# Creating Multivariate Rings

Referring to the rings and fields chapters for more details, our focus in this section is on examples
which illustrate how to create multivariate polynomial rings and their elements. Of particular
interest to us is the introduction of decorated ring types which allow one to implement multivariate
polynomial rings with gradings and/or filtrations.

## Types

OSCAR provides types for dense univariate and sparse multivariate polynomials. The univariate
ring types belong to the abstract type `PolyRing{T}`, their elements have abstract type
`PolyRingElem{T}`. The multivariate ring types belong to the abstract type `MPolyRing{T}`,
their elements have abstract type `MPolyRingElem{T}`. Here, `T` is the element type
of the coefficient ring of the polynomial ring.

Multivariate rings with gradings and/or filtrations are modelled by objects of type
`MPolyRing_dec{T, S}  :< MPolyRing{T}`, with elements of type
`MPolyRingElem_dec{T, S}  :< MPolyRingElem{T}`. Here, `S` is the element type of the
multivariate ring, and  `T` is the element type of its coefficient ring as above.
	
## Constructors

The basic constructor below allows one to build multivariate polynomial rings:

```@julia
PolynomialRing(C::Ring, v::Vector{String}; ordering=:lex, cached = true)
```

Its return value is a tuple, say `R, vars`, consisting of a polynomial ring `R` with coefficient ring `C` and a vector `vars` of generators (variables) which print according to the strings in the vector `v` .
The input `ordering=:lex` refers to the lexicograpical monomial ordering which specifies the default way of storing and displaying polynomials in OSCAR  (terms are sorted in descending
order). The other possible choices are `:deglex` and `:degrevlex`. Gröbner bases, however, can be computed with respect to any monomial ordering. See the section on Gröbner bases.

!!! note
    Caching is used to ensure that a given ring constructed from given parameters is unique in the system. For example, there is only one ring of multivariate polynomials over  $\mathbb{Z}$ in the variables x, y, z with `ordering=:lex`.

###### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"])
typeof(R)
typeof(x)
S, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"])
R === S
```

```@repl oscar
R1, x = PolynomialRing(QQ, ["x"])
typeof(x)
R2, (x,) = PolynomialRing(QQ, ["x"])
typeof(x)
R3, x = PolynomialRing(QQ, "x")
typeof(x)
```

```@repl oscar
V = ["x[1]", "x[2]"]
T, x = PolynomialRing(GF(3), V)
x
```

The constructor illustrated below allows for the convenient handling of variables with multi-indices:

```@repl oscar
R, x, y, z = PolynomialRing(QQ, "x" => (1:3, 1:4), "y" => 1:2, "z" => (1:1, 1:1, 1:1))
x
y
z
```

## Coefficient Rings 

Gröbner bases are implemented for multivariate polynomial rings over the fields and rings from this list:

###### The field of rational numbers $\mathbb{Q}$

```@repl oscar
QQ
```
###### Finite fields $\mathbb{F_p}$, $p$ a prime

```@repl oscar
GF(3)
GF(ZZ(2)^127 - 1)
```

###### Finite fields $\mathbb{F}_{p^n}$ with $p^n$ elements, $p$ a prime

```@repl oscar
FiniteField(2, 70, "a")
```

###### Simple algebraic extensions of $\mathbb{Q}$ or $\mathbb{F}_p$
  
```@repl oscar
T, t = PolynomialRing(QQ, "t")
K, a = NumberField(t^2 + 1, "a")
F = GF(3)
T, t = PolynomialRing(F, "t")
K, a = FiniteField(t^2 + 1, "a")
```

###### Purely transcendental extensions of $\mathbb{Q}$ or $\mathbb{F}_p$

```@repl oscar
T, t = PolynomialRing(QQ, "t")
QT = FractionField(T)
parent(t)
parent(1//t)
T, (s, t) = PolynomialRing(GF(3), ["s", "t"]);
QT = FractionField(T)
```

###### The ring of integers $\mathbb{Z}$

```@repl oscar
ZZ
```

## Gradings and Filtrations

The following functions implement multivariate polynomial rings decorated with $\mathbb Z$-gradings by (weighted) degree:

```@docs
grade(R::MPolyRing, W::Vector{Int})
```

```@docs
GradedPolynomialRing(C::Ring, V::Vector{String}, W::Vector{Int}; ordering=:lex)
```

!!! note
    OSCAR functionality for working with gradings by arbitrary abelian groups and for working with filtrations is currently under development.

## Data Associated to Multivariate Rings

If `R` is  a multivariate polynomial ring with coefficient ring `C`, then

- `base_ring(R)` refers to `C`,
- `gens(R)` to the generators (variables) of `R`,
- `ngens(R)` to the number of these generators, and
- `gen(R, i)` as well as `R[i]` to the `i`-th generator.

###### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
base_ring(R)
gens(R)
gen(R, 2)
R[3] 
ngens(R)
```

## Elements of Multivariate Rings

One way to create elements of a multivariate  polynomial ring, `R`, is
to build up polynomials from the generators (variables) of `R` using
basic arithmetic as shown below:

###### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
f = 3*x^2+y*z
typeof(f)
```

```@repl oscar
R, (x, y, z) = GradedPolynomialRing(QQ, ["x", "y", "z"])
f = 3*x^2+y*z
typeof(f)
```

Alternatively, there is the following constructor:

```@julia
(R::MPolyRing{T})(c::Vector{T}, e::Vector{Vector{Int}}) where T <: RingElem
```

Its return value is the element of  `R`  whose nonzero coefficients are specified by the elements of `c`,
with exponent vectors given by the elements of `e`.

###### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
f = 3*x^2+y*z
g = R(QQ.([3, 1]), [[2, 0, 0], [0, 1, 1]])
f == g
```

!!! note
    An often more effective way to create polynomials is to use the `MPoly` build context, see the chapter on rings.

Given an element `f` of a multivariate polynomial ring `R`,
- `parent(f)` refers to `R`,
- `monomial(f, i)` to the `i`-th monomial of `f`, 
- `term(f, i)` to the `i`-th term of `f`,
- `coeff(f, i)` to the coefficient of the `i`-th term of `f`, and
- `exponent_vector(f, i)` to the exponent vector of the `i`-th term of `f`.


###### Examples

```@repl oscar
R, (x, y) = PolynomialRing(GF(5), ["x", "y"])
c = map(GF(5), [1, 2, 3])
e = [[3, 2], [1, 0], [0, 1]]
f = R(c, e)
parent(f)
coeff(f, 2)
exponent_vector(f, 2)
monomial(f, 2)
term(f, 2)
```

## Homomorphisms of Multivariate Rings

How to handle homomorphisms of multivariate polynomial rings and their decorated types is described in
a more general context in the section on affine algebras. 
