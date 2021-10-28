```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["rings.md"]
```

# Creating Polynomial Rings

Referring to the corrresponding section of the ring chapter for details, we 
recall how to create multivariate polynomial rings. Then we show how to assign
gradings to these rings and how to set up maps between them.

## Construction

The standard constructor below allows one to build multivariate polynomial rings:

```@julia
PolynomialRing(C::Ring, v::Vector{String}; ordering=:lex)
```

Its return value is a tuple, say `R, vars`, consisting of a polynomial ring `R` with coefficient ring `C` and a vector `vars` of generators (variables) which print according to the strings in the vector `v` .
The input `ordering=:lex` refers to the lexicograpical monomial ordering which specifies the default way of storing and displaying polynomials in OSCAR  (terms are sorted in descending
order). The other possible choices are `:deglex` and `:degrevlex`. Gröbner bases, however, can be computed with respect to any monomial ordering. See the section on Gröbner bases.

!!! note
    The abstract type of multivariate polynomial rings is `MPolyRing`, while that of
    univariate polynomial rings is `PolyRing`. The elements of these rings are stored
    in a sparse and dense format,  respectively. 


###### Examples


```@repl oscar
Qx, x = PolynomialRing(QQ, ["x"])
typeof(x)
Qx, (x,) = PolynomialRing(QQ, ["x"])
typeof(x)
Qx, x = PolynomialRing(QQ, "x")
typeof(x)
```

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
f = x^2+y^2*z
typeof(f)
```

```@repl oscar
S, (a,b) = PolynomialRing(ZZ, ["a", "b"])
g = a*b
typeof(g)
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
```

## Coefficient Rings 

Gröbner bases are implemented for multivariate polynomial rings over fields and rings from this list:

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

## Data Associated to Polynomial Rings

If `R` is  a multivariate polynomial ring `R` with coefficient ring `C`, then

- `base_ring(R)` refers to `C`,
- `gens(R)` to the generators (variables) of `R`,
- `ngens(Q)` to the number of these generators, and
- `gen(R, i)` as well as `R[i]` to the `i`th generator.

###### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
base_ring(R)
gens(R)
gen(R, 2)
R[3] 
ngens(R)
```

## Elements of Polynomial Rings

One way to construct elements of a multivariate  polynomial ring `R` is
to build up polynomials from the generators (variables) of `R` using
basic arithmetic as shown below:

###### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
f = 3*x^2+y*z
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
f = R(QQ.([3, 1]), [[2, 0, 0], [0, 1, 1]])
QQ.([3, 1]) == map(QQ, [3, 1]) 
```

If an element `f` of `R` has been constructed, then
- `parent(f)` refers to `R`. 


###### Examples

```@repl oscar
R, (x, y) = PolynomialRing(GF(5), ["x", "y"])
c = map(GF(5), [1, 2, 3])
e = [[3, 2], [1, 0], [0, 1]]
f = R(c, e)
parent(f)
```

!!! note
    An often more effective way to construct polynomials  is to use the `MPoly` build context as explained in the chapter on rings.


## Gradings

The following functions allow one to assign gradings to multivariate polynomial rings:

```@docs
grade(R::MPolyRing, W::Vector{Int})
```

More directly, graded polynomial rings can be constructed as follows:

```@docs
GradedPolynomialRing(C::Ring, V::Vector{String}, W::Vector{Int}; ordering=:lex)
```

!!! note
    The return types of the constructors above are all subtypes of `MPolyRing`.

## Homomorphisms of Polynomial Rings

Functionality for dealing with homomorphisms of multivariate polynomial rings is described in the more general context of affine algebras.
