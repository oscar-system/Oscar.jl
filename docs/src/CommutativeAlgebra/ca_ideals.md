```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["ca_ideals.md"]
```

# Ideals

## Constructors

```@docs
ideal(g::Array{T, 1}) where {T <: MPolyElem}
```

## Gröbner Bases

### Monomial Orderings

### Normal Forms

### Computing Gröbner Bases

```@docs
groebner_basis(I::MPolyIdeal)
```

```@docs
groebner_basis(I::MPolyIdeal, ord::Symbol; complete_reduction::Bool=false)
```

#### Gröbner Bases with transformation matrix

```@docs
groebner_basis_with_transformation_matrix(I::MPolyIdeal; ord::Symbol = :degrevlex, complete_reduction::Bool=false)
```

    fglm

    Gröbner walks

    Hilbert-driven

#### Leading Ideals

```@docs
leading_ideal(g::Array{T, 1}, args...) where { T <: MPolyElem }
```

```@docs
leading_ideal(I::MPolyIdeal)
```

#### Gröbner Bases over the integers

#### ....


### Syzygies

#### Generators of syzygies

```@docs
syzygy_generators(a::Array{<:MPolyElem, 1})
```

## Data Associated to Ideals

```@docs
base_ring(I::MPolyIdeal)
```

### Number of Generators

```@docs
ngens(I::MPolyIdeal)
```

### Generators

```@docs
gens(I::MPolyIdeal)
```

### Dimension

```@docs
dim(I::MPolyIdeal)
```

### Codimension

```@docs
codim(I::MPolyIdeal)
```
    
## Operations on Ideals

### Simple ideal Operations

#### Powers of Ideal

```@docs
:^(I::MPolyIdeal, m::Int)
```
#### Sum of Ideals

```@docs
:+(I::MPolyIdeal, J::MPolyIdeal)
```

#### Product of Ideals

```@docs
:*(I::MPolyIdeal, J::MPolyIdeal)
```

### Intersection of Ideals

```@docs
intersect(I::MPolyIdeal, J::MPolyIdeal)
```

### Ideal Quotients

Given two ideals $I, J$ of a ring $R$, the ideal quotient of $I$ by $J$ is the ideal

$I:J= \bigl\{f \in R\:\big|\: f J \subset I\bigr\}\subset R.$

```@docs
quotient(I::MPolyIdeal, J::MPolyIdeal)
```

### Saturation

Given two ideals $I, J$ of a ring $R$, the saturation of $I$ with respect to $J$ is the ideal

$I:J^{\infty} = \bigl\{ f \in R \:\big|\: f J^k \!\subset I {\text{ for some }}k\geq 1 \bigr\} = \textstyle{\bigcup\limits_{k=1}^{\infty} (I:J^k)}.$

```@docs
saturation(I::MPolyIdeal, J::MPolyIdeal)
```

```@docs
saturation_with_index(I::MPolyIdeal, J::MPolyIdeal)
```

### Elimination

```@docs
eliminate(I::MPolyIdeal, lv::Array{MPolyElem, 1})
```

### Homogenization and Dehomogenization

    homogenize(I,t)   CAVEAT: Als Ideal! Auch poly, vector, etc. Siehe M2.

    dehomogenize(I,t)

## Tests on Ideals

### Equality of Ideals

```@docs
:(==)(I::MPolyIdeal, J::MPolyIdeal)
```

### Containment of Ideals

```@docs
issubset(I::MPolyIdeal, J::MPolyIdeal)
```

### Ideal Membership

```@docs
ideal_membership(f::MPolyElem, I::MPolyIdeal)
```

### Radical Membership

```@docs
radical_membership(f::MPolyElem, I::MPolyIdeal)
```

### Primality Test

```@docs
isprime(I::MPolyIdeal)
```

### Primary Test

```@docs
isprimary(I::MPolyIdeal)
```
  
    homogeneity test:        ishomogeneous(I)
    
    .....
	
## Decomposition of Ideals

### Radical

```@docs
radical(I::MPolyIdeal)
```

### Primary Decomposition

```@docs
primary_decomposition(I::MPolyIdeal)
```

### Minimal Associated Primes

```@docs
minimal_primes(I::MPolyIdeal)
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

### Absolute Primary Decomposition

    absolute_primary_decomposition(I)               --->  absPrimdecGTZ   

    


