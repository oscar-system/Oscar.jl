```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

# Ideals in Multivariate Rings

## Types

The OSCAR type for ideals in multivariate polynomial rings is of parametrized form
`MPolyIdeal{T}`, where `T` is the element type of the polynomial ring.

## Constructors

```@docs
ideal(R::MPolyRing, g::Vector)
```

## Data Associated to Ideals

### Basic Data

If `I` is an ideal of a multivariate polynomial ring  `R`, then

- `base_ring(I)` refers to `R`,
- `gens(I)` to the generators of `I`,
- `ngens(I)` to the number of these generators, and
- `gen(I, k)` as well as `I[k]` to the `k`-th such generator.

###### Examples

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> I = ideal(R, [x, y])^2
ideal(x^2, x*y, y^2)

julia> base_ring(I)
Multivariate polynomial ring in 2 variables x, y
  over rational field

julia> gens(I)
3-element Vector{QQMPolyRingElem}:
 x^2
 x*y
 y^2

julia> ngens(I)
3

julia> gen(I, 2)
x*y

```

### Dimension

```@docs
dim(I::MPolyIdeal)
```

### Codimension

```@docs
codim(I::MPolyIdeal)
```

### Degree

```@docs
degree(I::MPolyIdeal)
```

### Minimal Sets of Generators

In the graded case, we have:

```@docs
minimal_generating_set(I::MPolyIdeal{<:MPolyDecRingElem})
```
    
## Operations on Ideals

### Simple Ideal Operations

#### Powers of Ideal

```@docs
^(I::MPolyIdeal, m::Int)
```
#### Sum of Ideals

```@docs
+(I::MPolyIdeal{T}, J::MPolyIdeal{T}) where T
```

#### Product of Ideals

```@docs
*(I::MPolyIdeal{T}, J::MPolyIdeal{T}) where T
```

### Intersection of Ideals

```@docs
intersect(I::MPolyIdeal{T}, Js::MPolyIdeal{T}...) where T
intersect(V::Vector{MPolyIdeal{T}}) where T
```

### Ideal Quotients

Given two ideals $I, J$ of a ring $R$, the ideal quotient of $I$ by $J$ is the ideal

$I:J= \bigl\{f \in R\:\big|\: f J \subset I\bigr\}\subset R.$

```@docs
quotient(I::MPolyIdeal{T}, J::MPolyIdeal{T}) where T
```

### Saturation

Given two ideals $I, J$ of a ring $R$, the saturation of $I$ with respect to $J$ is the ideal

$I:J^{\infty} = \bigl\{ f \in R \:\big|\: f J^k \!\subset I {\text{ for some }}k\geq 1 \bigr\} = \textstyle{\bigcup\limits_{k=1}^{\infty} (I:J^k)}.$

```@docs
saturation(I::MPolyIdeal{T}, J::MPolyIdeal{T}) where T
saturation_with_index(I::MPolyIdeal{T}, J::MPolyIdeal{T}) where T
```

### Elimination

```@docs
eliminate(I::MPolyIdeal{T}, V::Vector{T}) where T <: MPolyRingElem
```

## Tests on Ideals

### Basic Tests

```@docs
is_zero(I::MPolyIdeal)
```

```@docs
is_one(I::MPolyIdeal)
```

```@docs
is_monomial(f::MPolyRingElem)
```

### Containment of Ideals

```@docs
is_subset(I::MPolyIdeal{T}, J::MPolyIdeal{T}) where T
```

### Equality of Ideals

```@docs
==(I::MPolyIdeal{T}, J::MPolyIdeal{T}) where T
```

### Ideal Membership

```@docs
ideal_membership(f::T, I::MPolyIdeal{T}) where T
```

### Radical Membership

```@docs
radical_membership(f::T, I::MPolyIdeal{T}) where T
```

### Primality Test

```@docs
is_prime(I::MPolyIdeal)
```

### Primary Test

```@docs
is_primary(I::MPolyIdeal)
```

## Decomposition of Ideals

We discuss various decomposition techniques. They are implemented for
polynomial rings over fields and, if explicitly mentioned, also for
polynomial rings over the integers. See [DGP99](@cite) for a survey.

### Radical

```@docs
radical(I::MPolyIdeal)
```

### Primary Decomposition

```@docs
primary_decomposition(I::MPolyIdeal; algorithm::Symbol = :GTZ)
```

### Absolute Primary Decomposition

```@docs
absolute_primary_decomposition(I::MPolyIdeal{QQMPolyRingElem})
```

### Minimal Associated Primes

```@docs
minimal_primes(I::MPolyIdeal; algorithm::Symbol = :GTZ)
```

### Weak Equidimensional Decomposition

```@docs
equidimensional_decomposition_weak(I::MPolyIdeal)
```

### Equidimensional Decomposition of radical

```@docs
equidimensional_decomposition_radical(I::MPolyIdeal)
```

### Equidimensional Hull

```@docs
equidimensional_hull(I::MPolyIdeal)
```

### Radical of the Equidimensional Hull

```@docs
equidimensional_hull_radical(I::MPolyIdeal)
```

## Homogenization and Dehomogenization

Referring to [KR05](@cite) for definitions and technical details, we discuss homogenization and dehomogenization in the context of $\mathbb Z^m$-gradings. 

```@docs
homogenization(f::MPolyRingElem, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName; pos::Union{Int,Nothing}=nothing)
```

```@docs
homogenization(f::MPolyRingElem, var::VarName; pos::Union{Int,Nothing}=nothing)
```

```@docs
dehomogenization(F::MPolyDecRingElem, pos::Int)
```


## Generating Special Ideals

### Katsura-n

These systems appeared in a problem of magnetism in physics.
For a given $n$ `katsura(n)` has $2^n$ solutions and is defined in a
polynomial ring with $n+1$ variables over the rational numbers. For a
given polynomial ring `R` with $n$ variables `katsura(R)` defines the
corresponding system with $2^{n-1}$ solutions.

```@docs
katsura(n::Int)
```

```@docs
katsura(R::MPolyRing)
```
