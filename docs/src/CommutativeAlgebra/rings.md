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

Its return value is a tuple, say `R, vars`, consisting of a polynomial ring `R` with coefficient ring `C` and an array `vars` of generators (variables) which print according to the strings in the supplied vector `v` .
The input `ordering=:lex` refers to the lexicograpical monomial ordering. This is the default monomial ordering in OSCAR for storing and printing polynomials. Other choices here  are `:deglex` and `:degrevlex`.
Gröbner bases, however, can be computed with respect to any monomial ordering. See the section on Gröbner bases.

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
f = x^2+y*z
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

If `f` is an element of `R`, then
- `parent(f)` refers to `R`. 


###### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
base_ring(R)
gens(R)
gen(R, 2)
R[3] 
ngens(R)
f = x^2+y*z
parent(f)
```


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
