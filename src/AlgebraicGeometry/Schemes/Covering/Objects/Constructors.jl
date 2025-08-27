########################################################################
# Constructors for Covering                                            #
########################################################################
### The default constructor
# Returns a scheme in which every affine patch is only
# glued to itself via the identity.
@doc raw"""
    Covering(patches::Vector{<:AbsAffineScheme})

Return a `Covering` with pairwise disjoint affine charts ``Uᵢ`` given by
the entries of `patches`. This `Covering` will have no gluings except
those gluings along the identity of every affine chart to itself.

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
```
"""
function Covering(patches::Vector{<:AbsAffineScheme})
  g = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsGluing}()
  for X in patches
    U = PrincipalOpenSubset(X)
    f = id_hom(U)
    g[X,X] = SimpleGluing(X, X, f, f, check=false)
  end
  return Covering(patches, g, check=false)
end

### Turns an affine scheme into a trivial covering
Covering(X::AbsAffineScheme) = Covering([X])

### The empty covering of the empty scheme over kk
empty_covering(kk::Ring) = Covering(kk)

@doc raw"""
    disjoint_union(C1::Covering, C2::Covering)

Return the `Covering` corresponding to the disjoint union of `C1` and `C2`.

The charts and gluings of the disjoint union are given by the disjoint union of the charts and gluings of the covers `C1` and `C2`.

# Examples
```jldoctest
julia> P1, (x,y) = QQ[:x, :y];

julia> P2, (u,v) = QQ[:u, :v];

julia> U1 = spec(P1);

julia> U2 = spec(P2);

julia> C1 = Covering(U1) # Set up the trivial covering with only one patch
Covering
  described by patches
    1: affine 2-space
  in the coordinate(s)
    1: [x, y]

julia> C2 = Covering(U2)
Covering
  described by patches
    1: affine 2-space
  in the coordinate(s)
    1: [u, v]

julia> C = disjoint_union(C1, C2)
Covering
  described by patches
    1: affine 2-space
    2: affine 2-space
  in the coordinate(s)
    1: [x, y]
    2: [u, v]
```
"""
function disjoint_union(C1::Covering, C2::Covering)
  C = Covering(vcat(patches(C1), patches(C2)))
  for (X, Y) in keys(gluings(C1))
    add_gluing!(C, C1[X, Y])
  end
  for (X, Y) in keys(gluings(C2))
    add_gluing!(C, C2[X, Y])
  end
  return C
end

