```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["quotient_rings.md"]
```

# Quotient Rings of Polynomial Rings

We discuss functionality for working with quotient rings of multivariate polynomial rings which may or may not have an assigned grading.

## Constructor

```@docs
quo(R::MPolyRing, I::MPolyIdeal)
```

!!! note
    With or without an assigned grading, the return types of `quo` are  subtypes of `MPolyQuo`.

!!! note
    In Oscar, elements of quotient rings are not necessarily reduced with regard to the modulus of the quotient ring.
    Operations involving Gröbner basis computations may lead to partial reductions. Full reductions, depending on the choice of a monomial ordering, are achieved by explicitly computing normal forms. The functions `simplify` and `simplify!` discussed in the sections below implements this.

## Data Associated to Quotient Rings of Polynomial Rings

### Basic Data

If `Q=R/I` is the quotient of a multivariate polynomial ring `R` modulo an ideal `I` of `R`, then

- `base_ring(Q)` refers to `R`,
- `modulus(Q)` to `I`,
- `gens(Q)` to the generators of `Q`, and
- `ngens(Q)` to the number of these generators.

###### Examples

```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
Q, _ = quo(R, ideal(R, [y-x^2, x-z^3]))
base_ring(Q)
modulus(Q)
gens(Q)
ngens(Q)
```

### Dimension

```@docs
dim(Q::MPolyQuo)
```

## Tests on Quotient Rings of Polynomial Rings

### Reducedness Test

```@docs
isreduced(Q::MPolyQuo)
```

## Elements of Quotient Rings

### Reducing Elements of Quotient Rings

```@docs
simplify(f::MPolyQuoElem)
```

### Tests on Elements of Quotient Rings

```@docs
 ==(f::MPolyQuoElem, g::MPolyQuoElem)
```

## Ideals in Quotient Rings of Polynomial Rings

### Constructors

```@docs
ideal(Q::MPolyQuo{T}, V::Vector{T}) where T <: MPolyElem
```

### Reducing Ideals in Quotient Rings

```@docs
simplify(a::MPolyQuoIdeal)
```

### Data Associated to Ideals in Quotient Rings

#### Basic Data

If `a` is an ideal of the quotient ring `Q`, then

- `base_ring(a)` refers to `Q`,
- `gens(a)` to the generators of `a`, and
- `ngens(a)` to the number of these generators.


#### Dimension of Ideals in Quotient Rings

```@docs
dim(a::MPolyQuoIdeal)
```

### Operations on Ideals in Quotient Rings 

#### Simple Ideal Operations in Quotient Rings

##### Powers of Ideal

```@docs
:^(a::MPolyQuoIdeal, m::Int)
```
##### Sum of Ideals

```@docs
:+(a::MPolyQuoIdeal, b::MPolyQuoIdeal)
```

##### Product of Ideals

```@docs
:*(a::MPolyQuoIdeal, b::MPolyQuoIdeal)
```

#### Intersection of Ideals

```@docs
intersect(a::MPolyQuoIdeal, bs::MPolyQuoIdeal...)
```

#### Ideal Quotients

```@docs
quotient(a::MPolyQuoIdeal, b::MPolyQuoIdeal)
```

### Tests on Ideals in Quotient Rings

#### Basic Tests

```@docs
iszero(a::MPolyQuoIdeal)
```

#### Equality of Ideals in Quotient Rings

```@docs
:(==)(a::MPolyQuoIdeal, b::MPolyQuoIdeal)
```

#### Containment of Ideals in Quotient Rings

```@docs
issubset(a::MPolyQuoIdeal, b::MPolyQuoIdeal)
```

