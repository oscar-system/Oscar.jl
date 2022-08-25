```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["creation.md"]
```

# Creating PBW-Algebras

## Types

PBW-algebras are modelled by objects of type `GAlgRing{T, S} <: NCRing`, their elements are objects of type
`GAlgElem{T, S} <: NCRingElem`. Here,  `T` is the element type of the field over which the PBW-algebra
is defined (the type `S` is added for internal use).

## Constructors

The basic constructor below allows one to build PBW-algebras:

```@docs
pbw_algebra(R::MPolyRing{T}, rel, ord::MonomialOrdering) where T
```

## Data Associated to PBW-Algebras

Given a PBW-algebra `A` over a field `K`, 

- `coefficient_ring(A)` refers to `K`,
- `gens(A)` to the generators of `A`,
- `ngens(A)` to the number of these generators, and
- `gen(A, i)` as well as `A[i]` to the `i`-th such generator.

###### Examples

```@repl oscar
R, (x,y,z) = QQ["x", "y", "z"];
L = [x*y, x*z, y*z + 1];
REL = strictly_upper_triangular_matrix(L);
A, (x,y,z) = pbw_algebra(R, REL, deglex(gens(R)));
coefficient_ring(A)
gens(A)
gen(A, 2)
A[3] 
ngens(A)
```

## Elements of PBW-Algebras

Elements of PBW-algebras are stored and printed in their standard representation.

### Constructors

One way to create elements of a PBW-algebra `A` over a field `K` is to build them up
from the generators of `A` using basic arithmetic as shown below:

###### Examples

```@repl oscar
R, (x,y,z) = QQ["x", "y", "z"];
L = [x*y, x*z, y*z + 1];
REL = strictly_upper_triangular_matrix(L);
A, (x,y,z) = pbw_algebra(R, REL, deglex(gens(R)));
f = 3*x^2+z*y
```

Alternatively, there is the following constructor:

```@julia
(A::PBWAlgRing)(cs::AbstractVector, es::AbstractVector{Vector{Int}})
```
Its return value is the element of  `A`  whose standard representation is built from
the elements of `cs` as coefficients, and the elements of `es` as exponents.

###### Examples

```@repl oscar
R, (x,y,z) = QQ["x", "y", "z"];
L = [x*y, x*z, y*z + 1];
REL = strictly_upper_triangular_matrix(L);
A, (x,y,z) = pbw_algebra(R, REL, deglex(gens(R)));
f = 3*x^2+z*y
g = A(QQ.([3, 1, 1]), [[2, 0, 0], [0, 1, 1], [0, 0, 0]])
f == g
```

An often more effective way to create elements is to use the corresponding build context as indicated below:

```@repl oscar
R, (x,y,z) = QQ["x", "y", "z"];
L = [x*y, x*z, y*z + 1];
REL = strictly_upper_triangular_matrix(L);
A, (x,y,z) = pbw_algebra(R, REL, deglex(gens(R)));
B =  build_ctx(A);
for i = 1:5 push_term!(B, QQ(i), [i+1, i, i-1]) end
finish(B)
```

### Special Elements

Given a PBW-algebra `A`, `zero(A)` and `one(A)` refer to the additive and multiplicative identity of `A`, respectively.
Relevant test calls on an element `f` of `A` are  `iszero(f)` and `isone(f)`.


### Data Associated to Elements of Multivariate Rings

Given an element `f` of a PBW-algebra `A`, 
- `parent(f)` refers to `A`,
