
########################################################################
# Constructors for CoveredScheme                                       #
########################################################################

### The standard constructor
@doc raw"""
    CoveredScheme(C::Covering)

Return a `CoveredScheme` ``X`` with `C` as its `default_covering`.


# Examples
```jldoctest
julia> P1, (x,y) = QQ["x", "y"];

julia> P2, (u,v) = QQ["u", "v"];

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

### Conversion of an affine scheme into a covered scheme
CoveredScheme(X::AbsAffineScheme) = CoveredScheme(Covering(X))

### Construct the empty covered scheme over the ring R
function empty_covered_scheme(R::RT) where {RT<:AbstractAlgebra.Ring}
  return CoveredScheme(R)
end

