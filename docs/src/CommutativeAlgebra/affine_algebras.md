```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["affine_algebras.md"]
```

# Affine Algebras

With regard to notation, we use *affine algebra* as a synonym for *quotient ring of a multivariate polynomial ring modulo an ideal*.
More specifically, if $R$ is a multivariate polynomial ring with coefficient ring $C$, and $A=R/I$ is the quotient ring of $R$
modulo an ideal $I$ of $R$, we refer to $A$ as an affine algebra over $C$, or an affine $C$-algebra. In this section, we discuss
functionality for handling such algebras in OSCAR.

!!! note
    As for the entire chapter on commutative algebra, most of the functions discussed here rely on Gröbner basis techniques. They are implemented for affine algebras over fields (exact fields supported by OSCAR) and, if not indicated otherwise, for affine algebras over the integers.

!!! note
    In Oscar, elements of quotient rings are not necessarily reduced with regard to the modulus of the quotient ring.
    Operations involving Gröbner basis computations may lead to partial reductions. Full reductions, depending on the choice of a monomial ordering, are achieved by explicitly computing normal forms. The functions `simplify` and `simplify!` discussed in this section implements this.


## Types

The OSCAR type for quotient rings of  multivariate polynomial rings is of parametrized form `MPolyQuo{T}`,
with elements of type `MPolyQuoElem{T}`. Here, `T` is the element type of the polynomial ring.
    
## Constructors

```@docs
quo(R::MPolyRing, I::MPolyIdeal)
```

## Data Associated to Affine Algebras

### Basic Data

If `A=R/I` is the quotient ring of a multivariate polynomial ring `R` modulo an ideal `I` of `R`, then

- `base_ring(A)` refers to `R`,
- `modulus(A)` to `I`,
- `gens(A)` to the generators of `A`, and
- `ngens(A)` to the number of these generators.

###### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
A, _ = quo(R, ideal(R, [y-x^2, z-x^3]))
base_ring(A)
modulus(A)
gens(A)
ngens(A)
```

### Dimension

```@docs
dim(A::MPolyQuo)
```

## Elements of Affine Algebras

### Types

The OSCAR type for elements of quotient rings of  multivariate polynomial rings is of
parametrized form `MPolyQuo{T}`, where `T` is the element type of the polynomial ring.

### Creating Elements of Affine Algebras

Elements of an affine algebra $R/I$ are created as images of elements of $R$ under the projection map.

###### Examples

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"]);
A, p = quo(R, ideal(R, [x^3*y^2-y^3*x^2, x*y^4-x*y^2]))
f = p(x^3*y^2-y^3*x^2+x*y)
typeof(f)
```

### Reducing Elements of Affine Algebras

```@docs
simplify(f::MPolyQuoElem)
```

### Tests on Elements of Affine Algebras

```@docs
==(f::MPolyQuoElem{T}, g::MPolyQuoElem{T}) where T
```

## Ideals in Affine Algebras

### Constructors

```@docs
ideal(Q::MPolyQuo{T}, V::Vector{T}) where T <: MPolyElem
```

### Reducing Ideals in Affine Algebras

```@docs
simplify(a::MPolyQuoIdeal)
```

### Data Associated to Ideals in Affine Algebras

#### Basic Data

If `a` is an ideal of the affine algebra `A`, then

- `base_ring(a)` refers to `A`,
- `gens(a)` to the generators of `a`, and
- `ngens(a)` to the number of these generators.

###### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);
A, _ = quo(R, ideal(R, [y-x^2, z-x^3]));
a = ideal(A, [x-y])
base_ring(a)
gens(a)
ngens(a)
```


#### Dimension of Ideals in Affine Algebras

```@docs
dim(a::MPolyQuoIdeal)
```

### Operations on Ideals in Affine Algebras

#### Simple Ideal Operations in Affine Algebras

##### Powers of Ideal

```@docs
:^(a::MPolyQuoIdeal, m::Int)
```
##### Sum of Ideals

```@docs
:+(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
```

##### Product of Ideals

```@docs
:*(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
```

#### Intersection of Ideals

```@docs
intersect(a::MPolyQuoIdeal{T}, bs::MPolyQuoIdeal{T}...) where T
```

#### Ideal Quotients

```@docs
quotient(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
```

### Tests on Ideals in Affine Algebras

#### Basic Tests

```@docs
iszero(a::MPolyQuoIdeal)
```

#### Equality of Ideals in Affine Algebras

```@docs
==(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
```

#### Containment of Ideals in Affine Algebras

```@docs
issubset(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
```

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

The usual methods for maps are supported, such
as `domain` and `codomain`.

```@docs
preimage(F::AlgHom, I::U) where U <: Union{MPolyIdeal, MPolyQuoIdeal}
kernel(F::AlgHom)
```

###### Examples

```@repl oscar
D1, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"])
C1, (s,t) = GradedPolynomialRing(QQ, ["s", "t"])
V1 = [s^3, s^2*t, s*t^2, t^3]
para = hom(D1, C1, V1)
twistedCubic = kernel(para)
C2, p2 = quo(D1, twistedCubic)
D2, (a, b, c) = GradedPolynomialRing(QQ, ["a", "b", "c"])
V2 = [p2(w-y), p2(x), p2(z)]
proj = hom(D2, C2, V2)
nodalCubic = kernel(proj)
```

```@repl oscar
D3,y = PolynomialRing(QQ, "y" => 1:3)
C3, x = PolynomialRing(QQ, "x" => 1:3)
V3 = [x[1]*x[2], x[1]*x[3], x[2]*x[3]]
F3 = hom(D3, C3, V3)
sphere = ideal(C3, [x[1]^3 + x[2]^3  + x[3]^3 - 1])
steinerRomanSurface = preimage(F3, sphere)
```

### Tests on Homomorphisms of Affine Algebras


```@docs
isinjective(F::AlgHom)
issurjective(F::AlgHom)
isbijective(F::AlgHom)
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

## Subalgebras

### Subalgebra Membership

```@docs
subalgebra_membership(f::T, v::Vector{T}) where T <: Union{MPolyElem, MPolyQuoElem}
```

### Minimal Subalgebra Generators

```@docs
minimal_subalgebra_generators(V::Vector{T}) where T <: Union{MPolyElem, MPolyQuoElem}
```

## Noether Normalization

```@docs
noether_normalization(A::MPolyQuo)
```

###### Example

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);
A, _ = quo(R, ideal(R, [x*y, x*z]));
L = noether_normalization(A);
L[1]
L[2]
L[3]
```

## Normalization of Rings

```@docs
normalization(A::MPolyQuo)
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

## Tests on Affine Algebras

### Reducedness Test

```@docs
isreduced(Q::MPolyQuo)
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

## Hilbert Series and Hilbert Polynomial

Let $R=K[x_1, \dots x_n]$ be a polynomial ring in $n$ variables over a field $K$.
Assign positive integer weights $w_i$ to the variables $x_i$, and grade
$R=\bigoplus_{d\geq 0} R_d$ according to the corresponding weighted degree. Let $I$ be an
ideal of $R$ which is homogeneous with respect to this
grading. Then the affine $K$-algebra $A=R/I$ inherits the grading:
$A = \bigoplus_{d\geq 0} A_d$, where each graded piece $A_d$ is a finite dimensional
$K$-vector space. In this situation, the *Hilbert function* of $A$ is
the function

$H(A, \underline{\phantom{d}}): \N \to \N, d \mapsto \dim_K(d).$

The *Hilbert series* of $A$ is the generating function

$H_A(t)=\sum_{d\geq 0} H(A, d) t^d.$

It can be written as a rational function in $t$, say, with denominator

$(1-t^{w_1})\cdots (1-t^{w_n}).$ 

Now suppose that the weights on the variables are all 1. Then we also have the *Hilbert
polynomial* $P_A(t)\in\mathbb{Q}[t]$ which satisfies $H(M,d)=P_M(d)$ for all $d \gg 0$.
Furthermore, the *degree* of $A$ is defined as the dimension of $A$ over $K$ if this dimension
is finite, and as the integer $d$ such that the leading term of the
Hilbert polynomial has the form $d t^e/e!$, otherwise.

!!! warning
    Currently, all functions described below are only implemented in the case where the weights on the variables are all 1.

```@docs
hilbert_series(A::MPolyQuo)
hilbert_series_reduced(A::MPolyQuo)
hilbert_series_expanded(A::MPolyQuo, d::Int)
hilbert_function(A::MPolyQuo, d::Int)
hilbert_polynomial(A::MPolyQuo)
degree(A::MPolyQuo)
```

###### Examples

```@repl oscar
R, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"]);
A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));
hilbert_series(A)
hilbert_series_reduced(A)
hilbert_series_expanded(A, 7)
hilbert_function(A,7)
hilbert_polynomial(A)
degree(A)
```
