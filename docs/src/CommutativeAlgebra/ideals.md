```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["ideals.md"]
```

# Ideals in Multivariate Rings

## Types

The OSCAR type for ideals in multivariate polynomial rings is of parametrized form
`MPolyIdeal{T}`, where `T` is the element type of the polynomial ring.

## Constructors

```@docs
ideal(g::Vector{T}) where {T <: MPolyElem}
```

## Gröbner Bases

### Computing Gröbner Bases

```@docs
groebner_basis(I::MPolyIdeal; ordering::Symbol = :degrevlex, complete_reduction::Bool = false)
```
```@docs
std_basis(I::MPolyIdeal, o::MonomialOrdering)
```
See e.g. [GP08](@cite) for the theoretical background on Gröbner- and standard bases.

#### Gröbner Bases with transformation matrix

```@docs
groebner_basis_with_transformation_matrix(I::MPolyIdeal; ordering::Symbol = :degrevlex, complete_reduction::Bool=false)
```

    fglm

    Gröbner walks

    Hilbert-driven

!!! warning "Expert functions for Gröbner bases"
    The following functions are low-level implementations of various Gröbner
    basis algorithms with many adjustable arguments. Only use these
    functions directly if you know what you are doing.

```@docs
f4( I::MPolyIdeal; initial_hts::Int=17, nr_thrds::Int=1, max_nr_pairs::Int=0, la_option::Int=2, reduce_gb::Int=1, info_level::Int=0)
```

#### Gröbner Bases over the integers

Over the integers the coefficients of the polynomials 
are not invertible, thus their handling when computing
Gröbner bases and normal forms plays an important role. This is done when 
computing strong Gröbner bases which ensure the following property: 
For any element of an ideal its leading term is divisible by a leading term of an 
element of a corresponding strong Gröbner basis.

The textbook [AL94](@cite) provides details on theory and algorithms as well as references.

```@repl oscar
R, (x,y) = PolynomialRing(ZZ, ["x","y"])
I = ideal(R, [2x,3x,4y])
H = groebner_basis(I)
```

### Leading Ideals

```@docs
leading_ideal(g::Vector{T}; ordering::MonomialOrdering) where { T <: MPolyElem }
leading_ideal(I::MPolyIdeal; ordering::MonomialOrdering)
```

### Normal Forms

```@docs
normal_form(f::T, J::MPolyIdeal) where { T <: MPolyElem }
normal_form(A::Vector{T}, J::MPolyIdeal) where { T <: MPolyElem }
```

### Syzygies

```@docs
syzygy_generators(a::Vector{<:MPolyElem})
```

## Data Associated to Ideals

### Basic Data

If `I` is an ideal of a multivariate polynomial ring  `R`, then

- `base_ring(I)` refers to `R`,
- `gens(I)` to the generators of `I`,
- `ngens(I)` to the number of these generators, and
- `gen(I, k)` as well as `I[k]` to the `k`-th such generator.

###### Examples

```@repl oscar
R, (x, y) = PolynomialRing(QQ, ["x", "y"])
I = ideal(R, [x, y])^2
base_ring(I)
gens(I)
ngens(I)
gen(I, 2)
```

### Dimension

```@docs
dim(I::MPolyIdeal)
```

### Codimension

```@docs
codim(I::MPolyIdeal)
```
### Minimal Sets of Generators

In the graded case, we have:

```@docs
minimal_generating_set(I::MPolyIdeal{<:MPolyElem_dec})
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
eliminate(I::MPolyIdeal{T}, l::Vector{T}) where T <: MPolyElem
```

## Tests on Ideals

### Basic Tests

```@docs
iszero(I::MPolyIdeal)
```

```@docs
isone(I::MPolyIdeal)
```

```@docs
is_monomial(f::MPolyElem)
```

### Containment of Ideals

```@docs
issubset(I::MPolyIdeal{T}, J::MPolyIdeal{T}) where T
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
primary_decomposition(I::MPolyIdeal; alg = :GTZ)
```

### Absolute Primary Decomposition

```@docs
absolute_primary_decomposition(I::MPolyIdeal{fmpq_mpoly})
```

### Minimal Associated Primes

```@docs
minimal_primes(I::MPolyIdeal; alg = :GTZ)
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
homogenization(f::MPolyElem, W::Union{fmpz_mat, Matrix{<:IntegerUnion}}, var::String, pos::Int = 1)
```

```@docs
homogenization(f::MPolyElem, var::String, pos::Int=1)
```

```@docs
dehomogenization(F::MPolyElem_dec, pos::Int)
```


## Generating Special Ideals

```@docs
katsura(n::Int)
```

```@docs
katsura(R::MPolyRing)
```
