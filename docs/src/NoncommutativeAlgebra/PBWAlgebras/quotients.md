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
Pages = ["quotients.md"]
```

# GR-Algebras: Quotients of PBW-Algebras

In analogy to the affine algebras section in the commutative algebra chapter, we describe OSCAR
functionality for dealing with quotients of PBW-algebras modulo two-sided ideals.

!!! note
    Quotients of PBW-algebras modulo two-sided ideals are also known as *GR-algebras* (here, *GR*
    stands for *Gröbner-Ready*; see [Lev05](@cite)).

## Types

GR-algebras are modelled by objects of type `PBWAlgQuo{T, S} <: NCRing`, their elements are objects of type
`PBWAlgQuoElem{T, S} <: NCRingElem`. Here,  `T` is the element type of the field over which the GR-algebra
is defined (the type `S` is added for internal use).


## Constructors

```@docs
quo(A::PBWAlgRing, I::PBWAlgIdeal)
```

### Exterior Algebras

The *$n$-th exterior algebra over a field $K$* is the quotient of the PBW-algebra

$A=K \langle e_1,\dots, e_n \mid e_ie_j = - e_je_i \ \text { for }\ i\neq j\rangle$

modulo the two-sided ideal

$\langle e_1^2,\dots, e_n^2\rangle.$

```@docs
    exterior_algebra(K::Ring, xs::Union{AbstractVector{<:AbstractString}, 
                                    AbstractVector{Symbol}, AbstractVector{Char}})
```

## Data Associated to Affine GR-Algebras

### Basic Data

If `Q=A/I` is the quotient ring of a PBW-algebra `A` modulo a two-sided ideal `I` of `A`, then

- `base_ring(Q)` refers to `A`,
- `modulus(Q)` to `I`,
- `gens(Q)` to the generators of `Q`,
- `ngens(Q)` to the number of these generators, and
- `gen(Q, i)` as well as `Q[i]` to the `i`-th such generator.

###### Examples

```jldoctest
julia> R, (x, y, z) = QQ["x", "y", "z"];

julia> L = [-x*y, -x*z, -y*z];

julia> REL = strictly_upper_triangular_matrix(L);

julia> A, (x, y, z) = pbw_algebra(R, REL, deglex(gens(R)));

julia> I = two_sided_ideal(A, [x^2, y^2, z^2]);

julia> Q, q = quo(A, I);

julia> base_ring(Q)
PBW-algebra over Rational Field in x, y, z with relations y*x = -x*y, z*x = -x*z, z*y = -y*z

julia> modulus(Q)
two_sided_ideal(x^2, y^2, z^2)

julia> gens(Q)
3-element Vector{PBWAlgQuoElem{fmpq, Singular.n_Q}}:
 x
 y
 z

julia> ngens(Q)
3

julia> gen(Q, 2)
y
```

## Elements of GR-Algebras

### Types

The OSCAR type for elements of quotient rings of  multivariate polynomial rings PBW-algebras is of
parametrized form `PBWAlgQuoElem{T, S}`, where `T` is the element type  of the
field over which the GR-algebra is defined (the type `S` is added for internal use).

### Creating Elements of GR-Algebras

Elements of a GR-algebra $Q = A/I$ are created as images of elements of $A$ under the projection map
or by directly coercing elements of $A$ into $Q$. The function `simplify!` reduces a given element
with regard to the modulus $I$.

###### Examples

```jldoctest
julia> R, (x, y, z) = QQ["x", "y", "z"];

julia> L = [-x*y, -x*z, -y*z];

julia> REL = strictly_upper_triangular_matrix(L);

julia> A, (x, y, z) = pbw_algebra(R, REL, deglex(gens(R)));

julia> I = two_sided_ideal(A, [x^2, y^2, z^2]);

julia> Q, q = quo(A, I);

julia> f = q(y*x+z^2)
-x*y + z^2

julia> typeof(f)
PBWAlgQuoElem{fmpq, Singular.n_Q}

julia> simplify!(f);

julia> f
-x*y

julia> g = Q(y*x+x^2)
x^2 - x*y

julia> f == g
true
```

### Data associated to Elements of GR-Algebras

Given an element `f` of an affine  GR-algebra `Q`, 

- `parent(f)` refers to `Q`.

## Ideals in GR-Algebras
