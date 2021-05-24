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
to construct multivariate polynomial rings.

## Admissible Rings of Coefficients

Here is a list of rings of coefficients which allow for Gröbner basis computations :

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

###### Pure transcendental extensions of $\mathbb{Q}$ or $\mathbb{F}_p$

###### The ring of integers $\mathbb{ZZ}$

```@repl oscar
ZZ
```

## Constructing Polynomial Rings


The basic constructor for multivariate polynomial rings reads as follows:

```@julia
PolynomialRing(C::Ring, v::Vector{String}; ordering=:lex)
```

The use of this constructor and that of additional constructors which allow for indexed variables is indicated below.

###### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
f = x^2+y*z
typeof(f)
S, (a,b) = PolynomialRing(ZZ, ["a", "b"])
g = a*b
typeof(g)
```

```@repl oscar
Qx, x = PolynomialRing(QQ, 2, "x");
x
```

```@repl oscar
Qx, x = @PolynomialRing(QQ, x[1:3]);
x
Qx, x = @PolynomialRing(QQ, x[1:3, 1:3]);
x
Qx, x, y = @PolynomialRing(QQ, x[1:3], y[1:2]);
x, y
I = [1, 3, 5];
Qx, x = @PolynomialRing(QQ, x[I]);
x
```

## Grading Polynomial Rings

```@docs
    grade(R::MPolyRing, v::Array{Int, 1})
```

!!! note
    
    
