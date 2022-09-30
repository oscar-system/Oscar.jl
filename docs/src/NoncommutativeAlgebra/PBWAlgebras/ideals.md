```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["ideals.md"]
```

# Ideals in PBW-algebras

## Types

The OSCAR type for ideals in PBW-algebras is of parametrized form
`PBWAlgIdeal{D, T, S}`, where `D` encodes the direction left, right,
or two-sided, and `T` is the element type of the field over which
the PBW-algebra is defined (the type `S` is added for internal use).

## Constructors

```@docs
left_ideal(g::Vector{<:PBWAlgElem})
right_ideal(g::Vector{<:PBWAlgElem})
two_sided_ideal(g::Vector{<:PBWAlgElem})
```

## Gröbner bases

## Data Associated to Ideals

If `I` is an ideal of a PBW-algebra  `A`, then

- `base_ring(I)` refers to `A`,
- `gens(I)` to the generators of `I`,
- `ngens(I)` to the number of these generators, and
- `gen(I, k)` as well as `I[k]` to the `k`-th such generator.

###### Examples

```@repl oscar
D, (x, y, dx, dy) = weyl_algebra(QQ, ["x", "y"])
I = left_ideal(D, [x, dx])
base_ring(I)
gens(I)
ngens(I)
gen(I, 2)
```

## Operations on Ideals

### Simple Ideal Operations

#### Powers of Ideal

#### Sum of Ideals

```@docs
+(I::PBWAlgIdeal{D, T, S}, J::PBWAlgIdeal{D, T, S}) where {D, T, S}
```

#### Product of Ideals

### Intersection of Ideals

```@docs
intersect(I::PBWAlgIdeal{D, T, S}, Js::PBWAlgIdeal{D, T, S}...) where {D, T, S}
```

### Elimination

## Tests on Ideals

```@docs
iszero(I:: PBWAlgIdeal)
```

```@docs
isone(I:: PBWAlgIdeal)
```

```@docs
issubset(I::PBWAlgIdeal{D, T, S}, J::PBWAlgIdeal{D, T, S}) where {D, T, S}
```

```@docs
==(I::PBWAlgIdeal{D, T, S}, J::PBWAlgIdeal{D, T, S}) where {D, T, S}
```

```@docs
ideal_membership(f::PBWAlgElem{T, S}, I::PBWAlgIdeal{D, T, S}) where {D, T, S}
```

