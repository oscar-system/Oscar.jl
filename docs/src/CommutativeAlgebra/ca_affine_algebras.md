```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["ca_affine_algebras.md"]
```

# Affine Algebras

An affine algebra over a field $K$, or an affine $K$-algebra, or simply an affine algebra, is the quotient $A=R/I$ of a multivariate
polynomial ring $R$ over $K$ modulo an ideal $I$ of $R$. The general constructor for quotients of multivariate polynomial rings
modulo ideals allows one to create such algebras.

## Constructor

```@docs
quo(R::MPolyRing, I::MPolyIdeal)
```

!!! note
    The return type of `quo` is a subtype of `MPolyQuo`.

## Data Associated to Affine Algebras

### Basic Data

If `A=R/I` is the quotient of a multivariate polynomial ring `R` modulo an ideal `I` of `R`, then

- `base_ring(A)` refers to `R`,
- `modulus(A)` to `I`,
- `gens(A)` to the generators of `A`, and
- `ngens(A)` to the number of these generators.

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
A, _ = quo(R, ideal(R, [y-x^2, x-z^3]))
base_ring(A)
modulus(A)
gens(A)
ngens(A)
```
i
### Dimension

```@docs
dim(A::MPolyQuo)
```

## Hilbert Series and Hilbert Polynomial

Let $R=K[x_1, \dots x_n]$ be a polynomial ring in $n$ variables over a field $K$.
Assign positive integer weights $w_i$ to the variables $x_i$, and grade
$R=\bigoplus_{d\geq 0} R_d$ according to the corresponding weighted degree. Let $I$ be an
ideal of $R$ which is homogeneous with respect to this
grading. Then the affine $K$-algebra $A=R/I$ inherits the grading:
$A = \bigoplus_{d\geq 0} A_d$, where each graded piece $A_d$ is a finite dimensional
$K$-vector space. In this situation, the *Hilbert function* of $A$ is
the function

$H(A, \underline{\phantom{d}}): \N \rightarrow \N, d \mapsto \dim_K(d).$

The *Hilbert series* of $A$ is the generating function

$H_A(t)=\sum_{d\geq 0} H(A, d) t^d.$

It can be written as a rational function in $t$, say, with denominator

$(1-t^{w_1})\cdots (1-t^{w_n}).$ 

Now suppose that the weights on the variables are all 1. Then we also have the *Hilbert
polynomial* $P_A(t)\in\mathbb{Q}[t]$ which satisfies $H(M,d)=P_M(d)$ for all $d \gg 0$.
Furthermore, the *degree* of $A$ is defined as the dimension of $A$ over $K$ if this dimension
is finite, and as the integer $d$ such that the leading term of the
Hilbert polynomial has the form $d t^e/e!$, otherwise.

CAVEAT: Currently only implemented in the case where the weights on the variables are all 1.

```@docs
hilbert_series(A::Union{MPolyRing, MPolyQuo})
```

```@docs
hilbert_series_reduced(A::Union{MPolyRing, MPolyQuo})
```

```@docs
hilbert_series_expanded(A::Union{MPolyRing, MPolyQuo}, d::Int)
```

```@docs
hilbert_function(A::Union{MPolyRing, MPolyQuo}, d::Int)
```

```@docs
hilbert_polynomial(A::Union{MPolyRing, MPolyQuo})
```

```@docs
degree(A::Union{MPolyRing, MPolyQuo})
```

###### Examples

```@repl oscar
R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
R,  (w, x, y, z) = grade(R)
A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));
hilbert_series(A)
hilbert_series_reduced(A)
hilbert_series_expanded(A, 7)
hilbert_function(A,7)
hilbert_polynomial(A)
degree(A)
```

## Ideals in Affine Algebras

### Constructors



## Tests on Affine Algebras

### Subalgebra Membership

```@docs
subalgebra_membership(f::S, v::Vector{S}) where S <: Union{MPolyElem{U}, MPolyQuoElem{U}} where U
```
```@repl oscar 
R, x = PolynomialRing(QQ, :x => 1:3)
f = x[1]^6*x[2]^6-x[1]^6*x[3]^6;
v = [x[1]^3*x[2]^3-x[1]^3*x[3]^3, x[1]^3*x[2]^3+x[1]^3*x[3]^3]
subalgebra_membership(f,v)
```

### Reducedness Test

```@docs
isreduced(A::MPolyQuo)
```

###### Example

```@repl oscar
R, (x,) = PolynomialRing(QQ, ["x"])
A, _ = quo(R, ideal(R, [x^4]))
isreduced(A)
```

### Normality Test

```@docs
isnormal(A::MPolyQuo)
```

###### Example

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
A, _ = quo(R, ideal(R, [z^2-x*y]))
isnormal(A)
```

#### Cohen-Macaulayness Test

iscohenmacaulay(R)

###### Example

## Homomorphisms of Affine Algebras

### Constructors

```@docs
AlgebraHomomorphism(D::U, C::W, V::Vector{X}) where 
{T, S <: MPolyElem{T},
U <: Union{MPolyRing{T}, MPolyQuo{S}},
W <: Union{MPolyRing{T}, MPolyQuo{S}},
X <: Union{S, MPolyQuoElem{S}}}
```

### Data Associated to Homomorphisms of Affine Algebras


```@docs
domain(F::AlgHom)
```

```@docs
codomain(F::AlgHom)
```

```@docs
preimage(F::AlgHom, I::U) where U <: Union{MPolyIdeal, MPolyQuoIdeal}
```

```@docs
kernel(F::AlgHom)
```

###### Examples

```@repl oscar
D1, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
C1, (s,t) = PolynomialRing(QQ, ["s", "t"])
V1 = [s^3, s^2*t, s*t^2, t^3]
para = hom(D1, C1, V1)
twistedCubic = kernel(para)
C2, _ = quo(D1, twistedCubic)
D2, (a, b, b) = PolynomialRing(QQ, ["a", "b", "c"])
V2 = [w-y, x, z]
proj = hom(D2, C2, V2)
nodalCubic = kernel(proj)
```

```@repl oscar
D3,y = PolynomialRing(QQ, :y => 1:3)
C3, x = PolynomialRing(QQ, :x => 1:3)
V3 = [x[1]*x[2], x[1]*x[3], x[2]*x[3]]
F3 = hom(D3, C3, V3)
sphere = ideal(C3, [x[1]^3 + x[2]^3  + x[3]^3 - 1])
steinerRomanSurface = preimage(F3, sphere)
```

### Tests on Homomorphisms of Affine Algebras

```@docs
isinjective(F::AlgHom)
```

```@docs
issurjective(F::AlgHom)
```

```@docs
isbijective(F::AlgHom)
```

```@docs
isfinite(F::AlgHom)
```

###### Examples

```@repl oscar
D, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
S, (a, b, c) = PolynomialRing(QQ, ["a", "b", "c"])
C, p = quo(S, ideal(S, [c-b^3]))
V = [p(2*a + b^6), p(7*b - a^2), p(c^2)]
F = hom(D, C, V)
issurjective(F)
D1, _ = quo(D, kernel(F))
F1 = hom(D1, C, V)
isbijective(F1)
```
```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, [ "x", "y", "z"])
C, (s, t) = PolynomialRing(QQ, ["s", "t"])
V = [s*t, t, s^2]
paraWhitneyUmbrella = hom(R, C, V)
D, _ = quo(R, kernel(paraWhitneyUmbrella))
isfinite(hom(D, C, V))
```

### Composition of Homomorphisms of Affine Algebras

```@docs
compose(F::AlgHom{T}, G::AlgHom{T}) where T
```

## Noether Normalization

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

## Normalization of Rings

```@docs
normalization(A::MPolyQuo)
```

```@docs
normalization_with_delta(A::MPolyQuo)
```

###### Examples

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
A, _ = quo(R, ideal(R, [(x^2-y^3)*(x^2+y^2)*x]))
L = normalization(A)
```

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
A, _ = quo(R, ideal(R, [z^3-x*y^4]))
L = normalization(A)
```

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
A, _ = quo(R, ideal(R, [(x^2-y^3)*(x^2+y^2)*x]))
L = normalization_with_delta(A)
```

## Integral Bases

```@docs
integral_basis(f::MPolyElem, i::Int)
```

###### Example


```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
f = (y^2-2)^2 + x^5
integral_basis(f, 2)
```
