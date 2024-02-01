```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

# Creating PBW-Algebras

## Types

PBW-algebras are modelled by objects of type `PBWAlgRing{T, S} <: NCRing`, their elements are objects of type
`PBWAlgElem{T, S} <: NCRingElem`. Here,  `T` is the element type of the field over which the PBW-algebra
is defined (the type `S` is added for internal use).

## Constructors

The basic constructor below allows one to build PBW-algebras:

```@docs
pbw_algebra(R::MPolyRing{T}, rel, ord::MonomialOrdering) where T
```
Some PBW-algebras are predefined in OSCAR.

### Weyl Algebras

The *$n$-th Weyl algebra over a field $K$* is the PBW-algebra
```math
D_n(K)=K \langle x_1,\ldots, x_n, \partial _1,\dots \partial _n \mid \partial_i x_i=x_i\partial _i +1, \partial _i x_j=x_j \partial _i \ \text { for }\ i\neq j\rangle.
```
Here,  we tacitly assume that
```math
x_j x_i=x_i x _j \; \text{ and } \; \partial _j \partial_i=\partial_i \partial _j \; \text{ for all } \; i,j.
```
Note that any  global monomial ordering on $\text{Mon}_{2n}(x, \partial)$ is admissible for $D_n(K)$.

The constructor below returns the algebras equipped with `degrevlex`.

```@docs
    weyl_algebra(K::Ring, xs::AbstractVector{<:VarName})
```

### Universal Enveloping Algebras of Finite Dimensional Lie Algebras

Let $\mathfrak g$ be an $n$-dimensional Lie algebra over a field $K$, and let $x_1, \dots, x_n$ be a $K$-basis of $\mathfrak g$.
Consider $n$ indeterminates which are also denoted by $x_1, \dots, x_n$.  The *universal enveloping algebra of $\mathfrak g$*
is the PBW-algebra
```math
U(\mathfrak g)=K \langle x_1,\ldots, x_n \mid x_jx_i = x_ix_j+[x_j, x_i],  \ 1\leq i<j \leq n \rangle,
```
where $[x_j, x_i]$ corresponds to evaluating the Lie bracket $[x_j, x_i]_\mathfrak g$. That the standard monomials
in $U(\mathfrak g)$ indeed form a $K$-basis for $U(\mathfrak g)$ is the content of the
Poincar``\'{\text{e}}``-Birkhoff-Witt theorem (the names PBW-basis and PBW-algebra are derived from this fact).

Note that any degree compatible global monomial ordering on $\N^n$ is admissible for $U(\mathfrak g)$.

The constructors below return the algebras equipped with `degrevlex`.

### Quantized Enveloping Algebras

### Non-Standard Quantum Deformation of $so_3$

## Data Associated to PBW-Algebras

Given a PBW-algebra `A` over a field `K`, 

- `coefficient_ring(A)` refers to `K`,
- `gens(A)` to the generators of `A`,
- `number_of_generators(A)` / `ngens(A)` to the number of these generators, and
- `gen(A, i)` as well as `A[i]` to the `i`-th such generator.

###### Examples

```jldoctest
julia> R, (x,y,z) = QQ["x", "y", "z"];

julia> L = [x*y, x*z, y*z + 1];

julia> REL = strictly_upper_triangular_matrix(L);

julia> A, (x,y,z) = pbw_algebra(R, REL, deglex(gens(R)));

julia> coefficient_ring(A)
Rational field

julia> gens(A)
3-element Vector{PBWAlgElem{QQFieldElem, Singular.n_Q}}:
 x
 y
 z

julia> gen(A, 2)
y

julia> A[3]
z 

julia> number_of_generators(A)
3

```

## Elements of PBW-Algebras

Elements of PBW-algebras are stored and printed in their standard representation.

### Constructors

One way to create elements of a PBW-algebra `A` over a field `K` is to build them up
from the generators of `A` using basic arithmetic as shown below:

###### Examples

```jldoctest
julia> R, (x,y,z) = QQ["x", "y", "z"];

julia> L = [x*y, x*z, y*z + 1];

julia> REL = strictly_upper_triangular_matrix(L);

julia> A, (x,y,z) = pbw_algebra(R, REL, deglex(gens(R)));

julia> f = 3*x^2+z*y
3*x^2 + y*z + 1

```

Alternatively, there is the following constructor:

```@julia
(A::PBWAlgRing)(cs::AbstractVector, es::AbstractVector{Vector{Int}})
```
Its return value is the element of  `A`  whose standard representation is built from
the elements of `cs` as coefficients, and the elements of `es` as exponents.

###### Examples

```jldoctest
julia> R, (x,y,z) = QQ["x", "y", "z"];

julia> L = [x*y, x*z, y*z + 1];

julia> REL = strictly_upper_triangular_matrix(L);

julia> A, (x,y,z) = pbw_algebra(R, REL, deglex(gens(R)));

julia> f = 3*x^2+z*y
3*x^2 + y*z + 1

julia> g = A(QQ.([3, 1, 1]), [[2, 0, 0], [0, 1, 1], [0, 0, 0]])
3*x^2 + y*z + 1

julia> f == g
true

```

An often more effective way to create elements is to use the corresponding build context as indicated below:

```jldoctest
julia> R, (x,y,z) = QQ["x", "y", "z"];

julia> L = [x*y, x*z, y*z + 1];

julia> REL = strictly_upper_triangular_matrix(L);

julia> A, (x,y,z) = pbw_algebra(R, REL, deglex(gens(R)));

julia> B =  build_ctx(A);

julia> for i = 1:5 push_term!(B, QQ(i), [i+1, i, i-1]) end

julia> finish(B)
5*x^6*y^5*z^4 + 4*x^5*y^4*z^3 + 3*x^4*y^3*z^2 + 2*x^3*y^2*z + x^2*y

```

### Special Elements

Given a PBW-algebra `A`, `zero(A)` and `one(A)` refer to the additive and multiplicative identity of `A`, respectively.
Relevant test calls on an element `f` of `A` are  `iszero(f)` and `isone(f)`.


### Data Associated to Elements of PBW-algebras

Given an element `f` of a PBW-algebra `A`, 
- `parent(f)` refers to `A`,

## Opposite Algebras

```@docs
opposite_algebra(A::PBWAlgRing)
```

If a map `opp` from a PBW-algebra to its opposite algebra is given,
then `inv(opp)` refers to the inverse of `opp`.
