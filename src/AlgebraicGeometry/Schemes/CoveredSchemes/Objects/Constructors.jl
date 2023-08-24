
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

julia> U1 = Spec(P1);

julia> U2 = Spec(P2);

julia> C = Covering([U1, U2]) # A Covering with two disjoint affine charts
Covering
  described by patches
    1: spec of multivariate polynomial ring
    2: spec of multivariate polynomial ring
  in the coordinate(s)
    1: [x, y]
    2: [u, v]

julia> V1 = PrincipalOpenSubset(U1, x); # Preparations for glueing

julia> V2 = PrincipalOpenSubset(U2, u);

julia> f = SpecMor(V1, V2, [1//x, y//x]); # The glueing isomorphism

julia> g = SpecMor(V2, V1, [1//u, v//u]); # and its inverse

julia> G = Glueing(U1, U2, f, g); # Construct the glueing

julia> add_glueing!(C, G) # Make the glueing part of the Covering
Covering
  described by patches
    1: spec of multivariate polynomial ring
    2: spec of multivariate polynomial ring
  in the coordinate(s)
    1: [x, y]
    2: [u, v]

julia> X = CoveredScheme(C) # Create a CoveredScheme from the Glueing
Scheme
  over rational field
with default covering
  described by patches
    1: spec of multivariate polynomial ring
    2: spec of multivariate polynomial ring
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
CoveredScheme(X::AbsSpec) = CoveredScheme(Covering(X))

### Construct the empty covered scheme over the ring R
function empty_covered_scheme(R::RT) where {RT<:AbstractAlgebra.Ring}
  return CoveredScheme(R)
end

