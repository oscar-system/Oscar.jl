```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["affine_algebras.md"]
```

# Affine Algebras and Their Ideals

With regard to notation, we use *affine algebra* as a synonym for *quotient ring of a multivariate polynomial ring modulo an ideal*.
More specifically, if $R$ is a multivariate polynomial ring with coefficient ring $C$, and $A=R/I$ is the quotient ring of $R$
modulo an ideal $I$ of $R$, we refer to $A$ as an *affine algebra over $C$*, or an *affine $C$-algebra*. In this section, we discuss
functionality for handling such algebras in OSCAR.

!!! note
    Most functions discussed here rely on Gröbner basis techniques. They are implemented for affine algebras over fields (exact fields supported by OSCAR) and, if not indicated otherwise, for affine algebras over the integers.

!!! note
    In OSCAR, elements of quotient rings are not necessarily reduced with regard to the modulus of the quotient ring.
    Operations involving Gröbner basis computations may lead to partial reductions. Full reductions, depending on the choice of a monomial ordering, are achieved by explicitly computing normal forms. The function `simplify` discussed in this section implements this.

!!! note
    Each grading on a multivariate polynomial ring `R`  in OSCAR  descends to a grading on the affine algebra `A = R/I`
    (recall that OSCAR ideals of graded polynomial rings are required to be homogeneous).
    Functionality for dealing with such gradings and our notation for describing this functionality descend accordingly.
	This applies, in particular, to the functions `ìs_graded`, `ìs_standard_graded`, `ìs_z_graded`, `ìs_zm_graded`,
	and `ìs_positively_graded` which will not be discussed again here. 

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
- `gens(A)` to the generators of `A`,
- `ngens(A)` to the number of these generators, and
- `gen(A, i)` as well as `A[i]` to the `i`-th such generator.

###### Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> A, _ = quo(R, ideal(R, [y-x^2, z-x^3]))
(Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z), Map from
Multivariate Polynomial Ring in x, y, z over Rational Field to Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z) defined by a julia-function with inverse)

julia> base_ring(A)
Multivariate Polynomial Ring in x, y, z over Rational Field

julia> modulus(A)
ideal(-x^2 + y, -x^3 + z)

julia> gens(A)
3-element Vector{MPolyQuoElem{fmpq_mpoly}}:
 x
 y
 z

julia> ngens(A)
3

julia> gen(A, 2)
y

```

In the graded case, we additionally have:

```@docs
grading_group(q::MPolyQuo{<:MPolyElem_dec})
```

### Dimension

```@docs
dim(A::MPolyQuo)
```

```@docs
vdim(A::MPolyQuo)
```

## Elements of Affine Algebras

### Types

The OSCAR type for elements of quotient rings of  multivariate polynomial rings is of
parametrized form `MPolyQuo{T}`, where `T` is the element type of the polynomial ring.

### Creating Elements of Affine Algebras

Elements of an affine algebra $A = R/I$ are created as images of elements of $R$ under the projection map
or by directly coercing elements of $R$ into $A$.

###### Examples

```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

julia> A, p = quo(R, ideal(R, [x^3*y^2-y^3*x^2, x*y^4-x*y^2]));

julia> f = p(x^3*y^2-y^3*x^2+x*y)
x^3*y^2 - x^2*y^3 + x*y

julia> typeof(f)
MPolyQuoElem{fmpq_mpoly}

julia> g = A(x^3*y^2-y^3*x^2+x*y)
x^3*y^2 - x^2*y^3 + x*y

julia> f == g
true

```

### Reducing Elements of Affine Algebras

```@docs
simplify(f::MPolyQuoElem)
```

### Tests on Elements of Affine Algebras

```@docs
==(f::MPolyQuoElem{T}, g::MPolyQuoElem{T}) where T
```

In the graded case, we additionally have:

```@docs
is_homogeneous(f::MPolyQuoElem{<:MPolyElem_dec})
```

### Data associated to Elements of Affine Algebras

Given an element `f` of an affine algebra `A`, 

- `parent(f)` refers to `A`.

In the graded case,  we also have:

```@docs
 homogeneous_components(f::MPolyQuoElem{<:MPolyElem_dec})
```

```@docs
homogeneous_component(f::MPolyQuoElem{<:MPolyElem_dec}, g::GrpAbFinGenElem)
```

```@docs
degree(f::MPolyQuoElem{<:MPolyElem_dec})
```

## Ideals in Affine Algebras

### Constructors

```@docs
ideal(Q::MPolyQuo{T}, V::Vector{T}) where T <: MPolyElem
```

### Reducing Generators of Ideals

```@docs
simplify(a::MPolyQuoIdeal)
```

### Data Associated to Ideals in Affine Algebras

#### Basic Data

If `a` is an ideal of the affine algebra `A`, then

- `base_ring(a)` refers to `A`,
- `gens(a)` to the generators of `a`,
- `ngens(a)` to the number of these generators,  and
- `gen(a, i)` as well as `a[i]` to the `i`-th such generator.

###### Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [y-x^2, z-x^3]));

julia> a = ideal(A, [x-y, z^4])
ideal(x - y, z^4)

julia> base_ring(a)
Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z)

julia> gens(a)
2-element Vector{MPolyQuoElem{fmpq_mpoly}}:
 x - y
 z^4

julia> ngens(a)
2

julia> gen(a, 2)
z^4

```

#### Dimension of Ideals in Affine Algebras

```@docs
dim(a::MPolyQuoIdeal)
```

#### Minimal Sets of Generators

In the graded case, we have:

```@docs
minimal_generating_set(I::MPolyQuoIdeal{<:MPolyElem_dec})
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

#### Containment of Ideals in Affine Algebras

```@docs
issubset(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
```

#### Equality of Ideals in Affine Algebras

```@docs
==(a::MPolyQuoIdeal{T}, b::MPolyQuoIdeal{T}) where T
```

#### Ideal Membership

```@docs
ideal_membership(f::MPolyQuoElem{T}, a::MPolyQuoIdeal{T}) where T
```

## Homomorphisms From Affine Algebras

If $A=R/I$ is an affine $C$-algebra, and $S$ is any ring, then defining a ring homomorphism
$\overline{\phi}: A \rightarrow S$ means to define a ring homomorphism $\phi: R \rightarrow S$
such that $I\subset \ker(\phi)$. Thus, $\overline{\phi} $ is determined by specifying its restriction
to $C$, and by assigning an image to each generator of $A$.
In OSCAR, such homomorphisms are created by using the following constructor:

```@docs
hom(A::MPolyQuo, S::NCRing, coeff_map, images::Vector; check::Bool = true)
```

Given a ring homomorphism `F` from `R` to `S` as above, `domain(F)` and `codomain(F)`
refer to `R` and `S`, respectively. Given ring homomorphisms `F` from `R` to `S` and
`G` from `S` to `T` as above, `compose(F, G)` refers to their composition.

## Homomorphisms of Affine Algebras

The OSCAR homomorphism type `AffAlgHom` models ring homomorphisms `R` $\to$ `S` such that
the type of both `R` and `S`  is a subtype of `Union{MPolyRing{T}, MPolyQuo{U}}`, where `T <: FieldElem` and
`U <: MPolyElem{T}`. Functionality for these homomorphism is discussed in what follows.
       
### Data Associated to Homomorphisms of Affine Algebras

```@docs
preimage(F::AffAlgHom, I::MPolyIdeal)
kernel(F::AffAlgHom)
```

###### Examples

```jldoctest
julia> D1, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"]);

julia> C1, (s,t) = GradedPolynomialRing(QQ, ["s", "t"]);

julia> V1 = [s^3, s^2*t, s*t^2, t^3];

julia> para = hom(D1, C1, V1)
Map with following data
Domain:
=======
Multivariate Polynomial Ring in w, x, y, z over Rational Field graded by
  w -> [1]
  x -> [1]
  y -> [1]
  z -> [1]
Codomain:
=========
Multivariate Polynomial Ring in s, t over Rational Field graded by
  s -> [1]
  t -> [1]

julia> twistedCubic = kernel(para)
ideal(-x*z + y^2, -w*z + x*y, -w*y + x^2)

julia> C2, p2 = quo(D1, twistedCubic);

julia> D2, (a, b, c) = GradedPolynomialRing(QQ, ["a", "b", "c"]);

julia> V2 = [p2(w-y), p2(x), p2(z)];

julia> proj = hom(D2, C2, V2)
Map with following data
Domain:
=======
Multivariate Polynomial Ring in a, b, c over Rational Field graded by
  a -> [1]
  b -> [1]
  c -> [1]
Codomain:
=========
Quotient of Multivariate Polynomial Ring in w, x, y, z over Rational Field graded by
  w -> [1]
  x -> [1]
  y -> [1]
  z -> [1] by ideal(-x*z + y^2, -w*z + x*y, -w*y + x^2)

julia> nodalCubic = kernel(proj)
ideal(-a^2*c + b^3 - 2*b^2*c + b*c^2)

```

```jldoctest
julia> D3,y = PolynomialRing(QQ, "y" => 1:3);

julia> C3, x = PolynomialRing(QQ, "x" => 1:3);

julia> V3 = [x[1]*x[2], x[1]*x[3], x[2]*x[3]];

julia> F3 = hom(D3, C3, V3)
Map with following data
Domain:
=======
Multivariate Polynomial Ring in y[1], y[2], y[3] over Rational Field
Codomain:
=========
Multivariate Polynomial Ring in x[1], x[2], x[3] over Rational Field

julia> sphere = ideal(C3, [x[1]^3 + x[2]^3  + x[3]^3 - 1])
ideal(x[1]^3 + x[2]^3 + x[3]^3 - 1)

julia> steinerRomanSurface = preimage(F3, sphere)
ideal(y[1]^6*y[2]^6 + 2*y[1]^6*y[2]^3*y[3]^3 + y[1]^6*y[3]^6 + 2*y[1]^3*y[2]^6*y[3]^3 + 2*y[1]^3*y[2]^3*y[3]^6 - y[1]^3*y[2]^3*y[3]^3 + y[2]^6*y[3]^6)

```

### Tests on Homomorphisms of Affine Algebras


```@docs
is_injective(F::AffAlgHom)
is_surjective(F::AffAlgHom)
is_bijective(F::AffAlgHom)
isfinite(F::AffAlgHom)
```

###### Examples

```jldoctest
julia> D, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> S, (a, b, c) = PolynomialRing(QQ, ["a", "b", "c"]);

julia> C, p = quo(S, ideal(S, [c-b^3]));

julia> V = [p(2*a + b^6), p(7*b - a^2), p(c^2)];

julia> F = hom(D, C, V)
Map with following data
Domain:
=======
Multivariate Polynomial Ring in x, y, z over Rational Field
Codomain:
=========
Quotient of Multivariate Polynomial Ring in a, b, c over Rational Field by ideal(-b^3 + c)

julia> is_surjective(F)
true

julia> D1, _ = quo(D, kernel(F));

julia> F1 = hom(D1, C, V);

julia> is_bijective(F1)
true

```

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, [ "x", "y", "z"]);

julia> C, (s, t) = PolynomialRing(QQ, ["s", "t"]);

julia> V = [s*t, t, s^2];

julia> paraWhitneyUmbrella = hom(R, C, V)
Map with following data
Domain:
=======
Multivariate Polynomial Ring in x, y, z over Rational Field
Codomain:
=========
Multivariate Polynomial Ring in s, t over Rational Field

julia> D, _ = quo(R, kernel(paraWhitneyUmbrella));

julia> isfinite(hom(D, C, V))
true

```

### Inverting Homomorphisms of Affine Algebras

```@docs
inverse(F::AffAlgHom)
```

```jldoctest
julia> D1, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> D, _ = quo(D1, [y-x^2, z-x^3])
(Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z), Map from
Multivariate Polynomial Ring in x, y, z over Rational Field to Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z) defined by a julia-function with inverse)

julia> C, (t,) = PolynomialRing(QQ, ["t"]);

julia> para = hom(D, C, [t, t^2, t^3]);

julia> is_bijective(para)
true

julia> inverse(para)
Map with following data
Domain:
=======
Multivariate Polynomial Ring in t over Rational Field
Codomain:
=========
Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z)

```

## Subalgebras

### Subalgebra Membership

```@docs
subalgebra_membership(f::T, V::Vector{T}) where T <: Union{MPolyElem, MPolyQuoElem}
```

### Minimal Subalgebra Generators

```@docs
minimal_subalgebra_generators(V::Vector{T}) where T <: Union{MPolyElem, MPolyQuoElem}
```

## Noether Normalization

```@docs
noether_normalization(A::MPolyQuo)
```

###### Examples

```jldoctest; setup = :(Singular.call_interpreter("""system("random", 47);"""))
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [x*y, x*z]));

julia> L = noether_normalization(A);

julia> L[1]
2-element Vector{MPolyQuoElem{fmpq_mpoly}}:
 -2*x + y
 -5*y + z

julia> L[2]
Map with following data
Domain:
=======
Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x*y, x*z)
Codomain:
=========
Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(2*x^2 + x*y, 10*x^2 + 5*x*y + x*z)

julia> L[3]
Map with following data
Domain:
=======
Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(2*x^2 + x*y, 10*x^2 + 5*x*y + x*z)
Codomain:
=========
Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x*y, x*z)

```
## Normalization

```@docs
normalization(A::MPolyQuo)
```

```@docs
normalization_with_delta(A::MPolyQuo)
```

## Integral Bases

```@docs
integral_basis(f::MPolyElem, i::Int)
```

## Tests on Affine Algebras

### Reducedness Test

```@docs
is_reduced(A::MPolyQuo)
```

### Normality Test

```@docs
is_normal(A::MPolyQuo)
```

### Cohen-Macaulayness Test

```@docs
is_cohen_macaulay(A::MPolyQuo)
```

## Hilbert Series and Hilbert Polynomial

Given a multivariate polynomial ring $R$ over a field $K$ together with a (multi)grading
on $R$ by a finitely generated abelian group $G$, let $I$ be an ideal of $R$ which is
homogeneous with respect to this grading. Then the affine $K-$algebra $A=R/I$
inherits the grading: $A = \bigoplus_{g\in G} A_g$. Suppose now that $R$ is positively
graded by $G$. That is, $G$ is free and each graded piece $R_g$ has finite dimension.
Then also $A_g$ is a finite dimensional $K$-vector space for each $g$, and we have the
well-defined *Hilbert function* of $A$,

$H(A, \underline{\phantom{d}}): G \to \N, \; g\mapsto \dim_K(A_g).$

The *Hilbert series* of $A$ is the generating function 

$H_A(\mathbb t)=\sum_{g\in G} H(A, g) \mathbb t^g$

(see  Section 8.2 in [MS05](@cite) for a formal discussion extending the classical case of
$\mathbb Z$-gradings with positive weights to the more general case of multigradings).
As in the classical case, the infinitely many values of the Hilbert function
can be expressed in finite terms by representing the Hilbert series as a rational function
(see Theorem 8.20 in [MS05](@cite) for a precise statement).

By a result of Macaulay, if $A = R/I$ is an affine algebra, and $L_{>}(I)$ is the leading
ideal of $I$ with respect to a global monomial ordering $>$, then the Hilbert function of $A$
equals that of $R/L_{>}(I)$ (see Theorem 15.26 in [Eis95](@cite)).
Thus, using Gröbner bases, the computation of Hilbert series can be reduced to the case where
the modulus of the affine algebra is a monomial ideal. In the latter case, we face a problem 
of combinatorial nature, and there are various strategies of how to proceed (see [KR05](@cite)).
The functions `hilbert_series`, `hilbert_series_reduced`, `hilbert_series_expanded`,
`hilbert_function`, `hilbert_polynomial`, and `degree` address the case of
$\mathbb Z$-gradings with positive weights, relying on corresponding Singular
functionality. The functions `multi_hilbert_series`, `multi_hilbert_series_reduced `,
and `multi_hilbert_function` offer a variety of different strategies and allow one to handle
positive gradings in general.

### $\mathbb Z$-Gradings With Positive Weights

Let $R=K[x_1, \dots x_n]$ be a polynomial ring in $n$ variables over a field $K$.
Assign positive integer weights $w_i$ to the variables $x_i$, and grade
$R=\bigoplus_{d\in \mathbb Z} R_d=\bigoplus_{d\geq 0} R_d$ according
to the corresponding weighted degree. Let $I$ be an
ideal of $R$ which is homogeneous with respect to this
grading. Then the affine $K$-algebra $A=R/I$ inherits the grading:
$A = \bigoplus_{d\geq 0} A_d$, where each graded piece $A_d$ is a finite dimensional
$K$-vector space. In this situation, the *Hilbert function* of $A$ is
of type

$H(A, \underline{\phantom{d}}): \N \to \N, \;d \mapsto \dim_K(d),$

and the *Hilbert series* of $A$ is the formal power series

$H_A(t)=\sum_{d\geq 0} H(A, d) t^d\in\mathbb Z[[t]].$

The Hilbert series can be written as a rational function $p(t)/q(t)$, with denominator

$q(t) = (1-t^{w_1})\cdots (1-t^{w_n}).$ 

In the standard $\mathbb Z$-graded case, where the weights on the variables are all 1, the Hilbert function is of polynomial nature: There exists
 a unique polynomial $P_A(t)\in\mathbb{Q}[t]$, the *Hilbert polynomial*, which satisfies $H(M,d)=P_M(d)$
for all $d \gg 0$. Furthermore, the *degree* of $A$ is defined as the dimension of $A$ over $K$ if this dimension
is finite, and as the integer $d$ such that the leading term of the Hilbert polynomial has the form $d t^e/e!$, otherwise.

```@docs
hilbert_series(A::MPolyQuo)
hilbert_series_reduced(A::MPolyQuo)
hilbert_series_expanded(A::MPolyQuo, d::Int)
hilbert_function(A::MPolyQuo, d::Int)
hilbert_polynomial(A::MPolyQuo)
degree(A::MPolyQuo)
```

### Positive Gradings in General

```@docs
multi_hilbert_series(A::MPolyQuo; alg::Symbol=:BayerStillmanA)
multi_hilbert_series_reduced(A::MPolyQuo; alg::Symbol=:BayerStillmanA)
multi_hilbert_function(A::MPolyQuo, g::GrpAbFinGenElem)
```
