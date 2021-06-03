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
to construct multivariate polynomial rings. Also, we show how to decorate such
rings with gradings.

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
T, t = PolynomialRing(QQ, :t=>1:2);
QT = FractionField(T)
parent(t[1])
parent(1//t[1])
U, u = PolynomialRing(GF(3), :u=>1:2);
QU = FractionField(U)
```

###### The ring of integers $\mathbb{ZZ}$

```@repl oscar
ZZ
```

## Constructions


The basic constructor for multivariate polynomial rings reads as follows:

```@julia
PolynomialRing(C::Ring, v::Vector{String}; ordering=:lex)
```

Its use and that of additional constructors which allow for indexed variables is indicated below.

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
v = ["x[1]", "x[2]"]
T, x = PolynomialRing(GF(3), v)
x
```

```@repl oscar
R, x = PolynomialRing(QQ, 2, "x");
x
x1
x[1]^2-x[2]
S, a, b, c= PolynomialRing(QQ, "a"=>1:3, "b"=>1:3, "c"=>1:5:10)
a, b, c
```

```@repl oscar
R, x = @PolynomialRing(QQ, x[1:3]);
x
x[1]^2-x[2]
R, x = @PolynomialRing(QQ, x[1:3, 1:3]);
x
R, x, y = @PolynomialRing(QQ, x[1:3], y[1:2]);
x, y
I = [1, 3, 5];
R, x = @PolynomialRing(QQ, x[I]);
x
```

## Gradings

```@docs
    grade(R::MPolyRing, v::Array{Int, 1})
```

    
    
