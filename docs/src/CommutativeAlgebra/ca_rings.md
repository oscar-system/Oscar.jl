```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["ca_rings.md"]
```

# Polynomial Rings

Referring to the sections on rings and fields for more details, we summarize how
to construct multivariate polynomial rings. Also, we show how to decorate
multivariate polynomial rings with gradings.

## Coefficients

Multivariate polynomial rings over coefficient rings from the list below allow for Gröbner basis computations :

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
Qx, x = QQ["x"]
K, a = NumberField(x^2 + 1, "a")
F = GF(3)
Fx, x = F["x"]
K, a = FiniteField(x^2 + 1, "a")
```

###### Purely transcendental extensions of $\mathbb{Q}$ or $\mathbb{F}_p$

```@repl oscar
T, (s, t) = PolynomialRing(QQ, ["s", "t"]);
QT = FractionField(T)
parent(t)
parent(1//t)
T, t = PolynomialRing(GF(3), "t");
QT = FractionField(T)
```

###### The ring of integers $\mathbb{Z}$

```@repl oscar
ZZ
```

## Constructions


The basic constructor for multivariate polynomial rings reads as follows:

```@julia
PolynomialRing(C::Ring, v::Vector{String}; ordering=:lex)
```

Its use, together with that of a convenient constructor for handling variables with multi-indices, is indicated below:

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

```@repl oscar
R, x, y, z = PolynomialRing(QQ, "x" => (1:3, 1:4), "y" => 1:2, "z" => (1:1, 1:1, 1:1))
x
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
    The return type of `PolynomialRing` as well as that of the constructors for graded rings is a subtype of `MPolyRing`.


