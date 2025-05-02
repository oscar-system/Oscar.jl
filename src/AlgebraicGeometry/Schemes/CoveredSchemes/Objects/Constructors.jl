
########################################################################
# Constructors for CoveredScheme                                       #
########################################################################

### The standard constructor
@doc raw"""
    CoveredScheme(C::Covering)

Return a `CoveredScheme` ``X`` with `C` as its `default_covering`.


# Examples
```jldoctest
julia> P1, (x,y) = QQ[:x, :y];

julia> P2, (u,v) = QQ[:u, :v];

julia> U1 = spec(P1);

julia> U2 = spec(P2);

julia> C = Covering([U1, U2]) # A Covering with two disjoint affine charts
Covering
  described by patches
    1: affine 2-space
    2: affine 2-space
  in the coordinate(s)
    1: [x, y]
    2: [u, v]

julia> V1 = PrincipalOpenSubset(U1, x); # Preparations for gluing

julia> V2 = PrincipalOpenSubset(U2, u);

julia> f = morphism(V1, V2, [1//x, y//x]); # The gluing isomorphism

julia> g = morphism(V2, V1, [1//u, v//u]); # and its inverse

julia> G = Gluing(U1, U2, f, g); # Construct the gluing

julia> add_gluing!(C, G) # Make the gluing part of the Covering
Covering
  described by patches
    1: affine 2-space
    2: affine 2-space
  in the coordinate(s)
    1: [x, y]
    2: [u, v]

julia> X = CoveredScheme(C) # Create a CoveredScheme from the Gluing
Scheme
  over rational field
with default covering
  described by patches
    1: affine 2-space
    2: affine 2-space
  in the coordinate(s)
    1: [x, y]
    2: [u, v]

```
"""
function CoveredScheme(C::Covering)
  refinements = Dict{Tuple{Covering, Covering}, CoveringMorphism}()
  X = CoveredScheme([C], refinements)
  set_attribute!(X, :seed_covering, C)
  return X
end

@doc raw"""
    disjoint_union(Xs::Vector{<:AbsCoveredScheme}) -> (AbsCoveredScheme, Vector{<:AbsCoveredSchemeMor})

Return the disjoint union of the non-empty vector of covered schemes as
a covered scheme.

# Input:
- a vector `Xs` of covered schemes.

# Output:
A pair ``(X, \mathrm{injections})`` where ``X`` is a covered scheme and
``\mathrm{injections}`` is a vector of inclusion morphisms ``ı_i\colon
X_i \to X``, where ``X`` is the disjoint union of the covered schemes
``X_i`` in `Xs`.

# Examples
```jldoctest
julia> R_1, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> I_1 = ideal(R_1, z*x^2 + y^3);

julia> X_1 = covered_scheme(proj(R_1, I_1))
Scheme
  over rational field
with default covering
  described by patches
    1: scheme((y//x)^3 + (z//x))
    2: scheme((x//y)^2*(z//y) + 1)
    3: scheme((x//z)^2 + (y//z)^3)
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> R_2, (u, v) = polynomial_ring(rational_field(), [:u, :v]);

julia> I_2 = ideal(R_2, u + v^2);

julia> X_2 = covered_scheme(spec(R_2, I_2))
Scheme
  over rational field
with default covering
  described by patches
    1: scheme(u + v^2)
  in the coordinate(s)
    1: [u, v]

julia> X, injections = disjoint_union([X_1, X_2]);

julia> X
Scheme
  over rational field
with default covering
  described by patches
    1: scheme((y//x)^3 + (z//x))
    2: scheme((x//y)^2*(z//y) + 1)
    3: scheme((x//z)^2 + (y//z)^3)
    4: scheme(u + v^2)
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]
    4: [u, v]

julia> injections
2-element Vector{CoveredSchemeMorphism{CoveredScheme{QQField}, CoveredScheme{QQField}, AbsAffineSchemeMor}}:
 Hom: scheme over QQ covered with 3 patches -> scheme over QQ covered with 4 patches
 Hom: scheme over QQ covered with 1 patch -> scheme over QQ covered with 4 patches
```
"""
function disjoint_union(Xs::Vector{<:AbsCoveredScheme})
  @req !is_empty(Xs) "Input should be a non-empty vector."
  covering = reduce(disjoint_union, default_covering.(Xs))
  X = CoveredScheme(covering)
  embed_covering(Xi) = Oscar.refinement_morphism(default_covering(Xi), covering)
  embed(Xi) = CoveredSchemeMorphism(Xi, X, embed_covering(Xi))
  injections = embed.(Xs)
  return X, injections
end

### Conversion of an affine scheme into a covered scheme
CoveredScheme(X::AbsAffineScheme) = CoveredScheme(Covering(X))

@doc raw"""
    covered_scheme(X::AbsAffineScheme) -> AbsCoveredScheme

Return a `CoveredScheme` ``C`` isomorphic to ``X``.

# Examples
```jldoctest
julia> R, x = polynomial_ring(rational_field(), :x => 1:3);

julia> X = spec(R);

julia> C = covered_scheme(X)
Scheme
  over rational field
with default covering
  described by patches
    1: affine 3-space
  in the coordinate(s)
    1: [x[1], x[2], x[3]]
```
"""
function covered_scheme(X::AbsAffineScheme)
  return CoveredScheme(X)
end

### Construct the empty covered scheme over the ring R
function empty_covered_scheme(R::RT) where {RT<:AbstractAlgebra.Ring}
  return CoveredScheme(R)
end
