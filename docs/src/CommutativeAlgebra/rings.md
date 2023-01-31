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
Pages = ["rings.md"]
```

# Creating Multivariate Rings

In this section, we illustrate by examples how to create multivariate polynomial rings and their elements,
while at the same time introducing and illustrating a special ring type for modelling multivariate polynomial
rings with gradings by finitely presented groups. For more details on multivariate polynomial rings, their
coefficient rings (fields), and their elements, we refer to the chapters on rings and fields. 

## Types

OSCAR provides types for dense univariate and sparse multivariate polynomials. The univariate
ring types belong to the abstract type `PolyRing{T}`, their elements have abstract type
`PolyRingElem{T}`. The multivariate ring types belong to the abstract type `MPolyRing{T}`,
their elements have abstract type `MPolyRingElem{T}`. Here, `T` is the element type
of the coefficient ring of the polynomial ring.

## Constructors

The basic constructor below allows one to build multivariate polynomial rings:

```@julia
PolynomialRing(C::Ring, V::Vector{String}; ordering=:lex, cached = true)
```

Its return value is a tuple, say `R, vars`, consisting of a polynomial ring `R` with coefficient ring `C` and a vector `vars` of generators (variables) which print according to the strings in the vector `V` .
The input `ordering=:lex` refers to the lexicograpical monomial ordering which specifies the default way of storing and displaying polynomials in OSCAR  (terms are sorted in descending
order). The other possible choices are `:deglex` and `:degrevlex`. Gröbner bases, however, can be computed with respect to any monomial ordering. See the section on Gröbner bases.

!!! note
    Caching is used to ensure that a given ring constructed from given parameters is unique in the system. For example, there is only one ring of multivariate polynomials over  $\mathbb{Z}$ in the variables x, y, z with `ordering=:lex`.

###### Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Integer Ring, fmpz_mpoly[x, y, z])

julia> typeof(R)
FmpzMPolyRing

julia> typeof(x)
fmpz_mpoly

julia> S, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Integer Ring, fmpz_mpoly[x, y, z])

julia> R === S
true

```

```jldoctest
julia> R1, x = PolynomialRing(QQ, ["x"])
(Multivariate Polynomial Ring in x over Rational Field, fmpq_mpoly[x])

julia> typeof(x)
Vector{fmpq_mpoly} (alias for Array{fmpq_mpoly, 1})

julia> R2, (x,) = PolynomialRing(QQ, ["x"])
(Multivariate Polynomial Ring in x over Rational Field, fmpq_mpoly[x])

julia> typeof(x)
fmpq_mpoly

julia> R3, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rational Field, x)

julia> typeof(x)
fmpq_poly

```

```jldoctest
julia> T, x = PolynomialRing(GF(3), ["x[1]", "x[2]"]);

julia> x
2-element Vector{gfp_mpoly}:
 x[1]
 x[2]

```

The constructor illustrated below allows for the convenient handling of variables with multi-indices:

```jldoctest
julia> R, x, y, z = PolynomialRing(QQ, "x" => (1:3, 1:4), "y" => 1:2, "z" => (1:1, 1:1, 1:1));

julia> x
3×4 Matrix{fmpq_mpoly}:
 x[1, 1]  x[1, 2]  x[1, 3]  x[1, 4]
 x[2, 1]  x[2, 2]  x[2, 3]  x[2, 4]
 x[3, 1]  x[3, 2]  x[3, 3]  x[3, 4]

julia> y
2-element Vector{fmpq_mpoly}:
 y[1]
 y[2]

julia> z
1×1×1 Array{fmpq_mpoly, 3}:
[:, :, 1] =
 z[1, 1, 1]

```

## Coefficient Rings 

Gröbner and standard bases are implemented for multivariate polynomial rings over the fields and rings below:

### The field of rational numbers $\mathbb{Q}$

```jldoctest
julia> QQ
Rational Field

```
### Finite fields $\mathbb{F_p}$, $p$ a prime

```jldoctest
julia> GF(3)
Galois field with characteristic 3

julia> GF(ZZ(2)^127 - 1)
Galois field with characteristic 170141183460469231731687303715884105727

```

### Finite fields $\mathbb{F}_{p^n}$ with $p^n$ elements, $p$ a prime

```jldoctest
julia> FiniteField(2, 70, "a")
(Finite field of degree 70 over F_2, a)

```

### Simple algebraic extensions of $\mathbb{Q}$ or $\mathbb{F}_p$
  
```jldoctest
julia> T, t = PolynomialRing(QQ, "t")
(Univariate Polynomial Ring in t over Rational Field, t)

julia> K, a = NumberField(t^2 + 1, "a")
(Number field over Rational Field with defining polynomial t^2 + 1, a)

julia> F = GF(3)
Galois field with characteristic 3

julia> T, t = PolynomialRing(F, "t")
(Univariate Polynomial Ring in t over Galois field with characteristic 3, t)

julia> K, a = FiniteField(t^2 + 1, "a")
(Finite field of degree 2 over F_3, a)

```

### Purely transcendental extensions of $\mathbb{Q}$ or $\mathbb{F}_p$

```jldoctest
julia> T, t = PolynomialRing(QQ, "t")
(Univariate Polynomial Ring in t over Rational Field, t)

julia> QT = FractionField(T)
Fraction field of Univariate Polynomial Ring in t over Rational Field

julia> parent(t)
Univariate Polynomial Ring in t over Rational Field

julia> parent(1//t)
Fraction field of Univariate Polynomial Ring in t over Rational Field

julia> T, (s, t) = PolynomialRing(GF(3), ["s", "t"]);

julia> QT = FractionField(T)
Fraction field of Multivariate Polynomial Ring in s, t over Galois field with characteristic 3

```

### The ring of integers $\mathbb{Z}$

```jldoctest
julia> ZZ
Integer Ring

```

## Gradings

Given a polynomial ring $R = C[x_1, \dots, x_n]$, we may endow $R$ with various gradings.
The *standard $\mathbb Z$-grading*  on $R$ is the decomposition
$R=\bigoplus_{d\in \mathbb Z} R_d=\bigoplus_{d\geq 0} R_d$ by the usual degree of polynomials.
More general $\mathbb Z$-gradings are obtained by assigning integer weights to the variables
and considering the corresponding weighted degrees. Even more generally, we may consider
multigradings: Given a finitely generated abelian group $G$, a *multigrading* on $R$ by $G$,
or a *$G$-grading*, or simply a *grading*, corresponds to a semigroup homomorphism
$\phi: \mathbb N^n \rightarrow G$: Given $\phi$, the *degree* of a monomial $x^\alpha$
is the image $\deg(x^\alpha):=\phi(\alpha)\in G$; the induced $G$-grading on $R$
is the decomposition $R = \bigoplus_{g\in G} R_g$ satisfying $R_g\cdot R_h\subset R_{g+h}$,
where $R_g$ is the free $C$-module generated by the monomials of degree $g$. This grading is determined by
assigning the *weights* $\deg(x_i)$ to the $x_i$. In other words, it is determined by  the *weight
vector* $W = (\deg(x_1), \dots, \deg(x_n))\in G^n.$ 

We refer to the textbooks [MS05](@cite) and [KR05](@cite) for details on multigradings. With respect to notation,
we follow the former book.

!!! note
    Given a $G$-grading on $R$, we refer to $G$ as the *grading group* of $R$. Moreover, we then say
    that $R$ is *$G$-graded*, or simply that $R$ is *graded*.
    If $R$ is a polynomial ring over a field, we say that a $G$-grading on $R$ is *positive* if $G$ is free
	and each graded part $R_g$, $g\in G$, has finite dimension. We then also say that $R$ is
	*positively graded (by $G$)*. Note that the positivity condition can be equivalently expressed by
	asking that $G$ is free and that the degree zero part consists of the constants only (see Theorem 8.6 in [MS05](@cite)).
		
!!! note
    Given a  `G`-grading on `R` in OSCAR, we say that `R` is *$\mathbb Z^m$-graded* if `is_free(G) && ngens(G) == rank(G) == m`
	evaluates to `true`. In this case, conversion routines allow one to switch back and forth between elements
	of `G` and integer vectors of length `m`. Specifically, if `R` is *$\mathbb Z$-graded*, that is,
	`is_free(G) && ngens(G) == rank(G) == 1` evaluates to `true`,  elements of `G` may be converted
	to integers and vice versa.

### Types

Multivariate rings with gradings are modelled by objects of type
`MPolyRing_dec{T, S}  :< MPolyRing{T}`, with elements of type
`MPolyRingElem_dec{T, S}  :< MPolyRingElem{T}`. Here, `S` is the element type of the
multivariate ring, and  `T` is the element type of its coefficient ring as above.

!!! note
    The types `MPolyRing_dec{T, S}` and `MPolyRingElem_dec{T, S}` are
    also meant to eventually model multivariate rings with filtrations
	and their elements.


The following function allows one to distinguish between graded and filtered rings:

```@docs
is_graded(R::MPolyRing_dec)
```

### Constructors for Graded Rings

There are two basic ways of creating multivariate rings with gradings:
While the `grade` function allows one to assign a grading to a polynomial ring already constructed,
the `GradedPolynomialRing` function is meant to create a graded polynomial ring all at once.

```@docs
grade(R::MPolyRing, W::Vector{GrpAbFinGenElem})
```

```@docs
grade(R::MPolyRing, W::Vector{<:Vector{<:IntegerUnion}})
```

```@docs
grade(R::MPolyRing, W::Vector{<:IntegerUnion})
```
```@docs
GradedPolynomialRing(C::Ring, V::Vector{String}, W; ordering=:lex)
```

## Tests on Graded Rings

```@docs
is_standard_graded(R::MPolyRing_dec)
```

```@docs
is_z_graded(R::MPolyRing_dec)
```

```@docs
is_zm_graded(R::MPolyRing_dec)
```

```@docs
is_positively_graded(R::MPolyRing_dec)
```

## Data Associated to Multivariate Rings

Given  a multivariate polynomial ring `R` with coefficient ring `C`, 

- `coefficient_ring(R)` refers to `C`,
- `gens(R)` to the generators (variables) of `R`,
- `ngens(R)` to the number of these generators, and
- `gen(R, i)` as well as `R[i]` to the `i`-th such generator.

###### Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> coefficient_ring(R)
Rational Field

julia> gens(R)
3-element Vector{fmpq_mpoly}:
 x
 y
 z

julia> gen(R, 2)
y

julia> R[3]
z 

julia> ngens(R)
3

```

In the graded case, we additionally have:

```@docs
grading_group(R::MPolyRing_dec)
```

```@docs
homogeneous_component(R::MPolyRing_dec, g::GrpAbFinGenElem)
```

## Elements of Multivariate Rings

### Constructors

One way to create elements of a multivariate  polynomial ring is
to build up polynomials from the generators (variables) of the ring using
basic arithmetic as shown below:

###### Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> f = 3*x^2+y*z
3*x^2 + y*z

julia> typeof(f)
fmpq_mpoly

julia> S, (x, y, z) = grade(R)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> g = 3*x^2+y*z
3*x^2 + y*z

julia> typeof(g)
MPolyElem_dec{fmpq, fmpq_mpoly}

julia> g == S(f)
true

```

Alternatively, there is the following constructor:

```@julia
(R::MPolyRing{T})(c::Vector{T}, e::Vector{Vector{Int}}) where T <: RingElem
```

Its return value is the element of  `R`  whose nonzero coefficients are specified by the elements of `c`,
with exponent vectors given by the elements of `e`.

###### Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> f = 3*x^2+y*z
3*x^2 + y*z

julia> g = R(QQ.([3, 1]), [[2, 0, 0], [0, 1, 1]])
3*x^2 + y*z

julia> f == g
true

```

An often more effective way to create polynomials is to use the `MPoly` build context as indicated below:

```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> B = MPolyBuildCtx(R)
Builder for an element of Multivariate Polynomial Ring in x, y over Rational Field

julia> for i = 1:5 push_term!(B, QQ(i), [i, i-1]) end

julia> finish(B)
5*x^5*y^4 + 4*x^4*y^3 + 3*x^3*y^2 + 2*x^2*y + x

```


### Special Elements

Given a multivariate polynomial ring `R`, `zero(R)` and `one(R)` refer to the additive and multiplicative identity of `R`, respectively.
Relevant test calls on an element `f` of `R` are  `iszero(f)` and `isone(f)`.


### Data Associated to Elements of Multivariate Rings

Given an element `f` of a multivariate polynomial ring `R` or a graded version of such a ring, 
- `parent(f)` refers to `R`,
- `total_degree(f)` to the total degree of `f`,
- `monomial(f, i)` to the `i`-th monomial of `f`, 
- `term(f, i)` to the `i`-th term of `f`,
- `coeff(f, i)` to the coefficient of the `i`-th term of `f`, and
- `exponent_vector(f, i)` to the exponent vector of the `i`-th term of `f`.


###### Examples

```jldoctest
julia> R, (x, y) = PolynomialRing(GF(5), ["x", "y"])
(Multivariate Polynomial Ring in x, y over Galois field with characteristic 5, gfp_mpoly[x, y])

julia> c = map(GF(5), [1, 2, 3])
3-element Vector{gfp_elem}:
 1
 2
 3

julia> e = [[3, 2], [1, 0], [0, 1]]
3-element Vector{Vector{Int64}}:
 [3, 2]
 [1, 0]
 [0, 1]

julia> f = R(c, e)
x^3*y^2 + 2*x + 3*y

julia> parent(f)
Multivariate Polynomial Ring in x, y over Galois field with characteristic 5

julia> total_degree(f)
5

julia> coeff(f, 2)
2

julia> exponent_vector(f, 2)
2-element Vector{Int64}:
 1
 0

julia> monomial(f, 2)
x

julia> term(f, 2)
2*x

```

Further functionality is available in the graded case:

```@docs
homogeneous_components(f::MPolyElem_dec{T, S}) where {T, S}
```

```@docs
homogeneous_component(f::MPolyElem_dec, g::GrpAbFinGenElem)
```

```@docs
is_homogeneous(f::MPolyElem_dec)
```

```@docs
degree(f::MPolyElem_dec)
```

## Homomorphisms From Multivariate Rings

If $R$ is a multivariate polynomial ring, and $S$ is any ring, then a ring homomorphism
$R \rightarrow S$ is determined by specifying its restriction to the coefficient ring of $R$,
and by assigning an image to each variable of $R$.
In OSCAR, such homomorphisms are created by using the following constructor:

```@docs
hom(R::MPolyRing, S::NCRing, coeff_map, images::Vector; check::Bool = true)
```

Given a ring homomorphism `F` from `R` to `S` as above, `domain(F)` and `codomain(F)`
refer to `R` and `S`, respectively.

!!! note
    The OSCAR homomorphism type `AffAlgHom` models ring homomorphisms `R` $\to$ `S` such that
    the type of both `R` and `S`  is a subtype of `Union{MPolyRing{T}, MPolyQuo{U}}`, where `T <: FieldElem` and
    `U <: MPolyElem{T}`. Functionality for these homomorphism is discussed in the section on affine algebras.

