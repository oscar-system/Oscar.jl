export empty_covering, disjoint_union, standard_covering
########################################################################
# Constructors for Covering                                            #
########################################################################
### The default constructor
# Returns a scheme in which every affine patch is only 
# glued to itself via the identity.
@Markdown.doc """
    Covering(patches::Vector{<:AbsSpec})

Return a `Covering` with pairwise disjoint affine charts ``Uᵢ`` given by 
the entries of `patches`. This `Covering` will have no glueings except 
those glueings along the identity of every affine chart to itself.

# Examples
```jldoctest
julia> P1, (x,y) = QQ["x", "y"];

julia> P2, (u,v) = QQ["u", "v"];

julia> U1 = Spec(P1);

julia> U2 = Spec(P2);

julia> C = Covering([U1, U2]) # A Covering with two disjoint affine charts
Covering with 2 patches
```
"""
function Covering(patches::Vector{<:AbsSpec})
  g = IdDict{Tuple{AbsSpec, AbsSpec}, AbsGlueing}()
  for X in patches
    U = PrincipalOpenSubset(X)
    f = identity_map(U)
    g[X,X] = SimpleGlueing(X, X, f, f, check=false)
  end
  return Covering(patches, g, check=false)
end

### Turns an affine scheme into a trivial covering
Covering(X::AbsSpec) = Covering([X])

### The empty covering of the empty scheme over kk
empty_covering(kk::Ring) = Covering(kk)

@Markdown.doc """
    disjoint_union(C1::Covering, C2::Covering)

Return the `Covering` corresponding to the disjoint union of `C1` and `C2`. 

The charts and glueings of the disjoint union are given by the disjoint union of the charts and glueings of the covers `C1` and `C2`.

# Examples
```jldoctest
julia> P1, (x,y) = QQ["x", "y"];

julia> P2, (u,v) = QQ["u", "v"];

julia> U1 = Spec(P1);

julia> U2 = Spec(P2);

julia> C1 = Covering(U1) # Set up the trivial covering with only one patch
Covering with 1 patch

julia> C2 = Covering(U2)
Covering with 1 patch

julia> C = disjoint_union(C1, C2)
Covering with 2 patches
```
"""
function disjoint_union(C1::Covering, C2::Covering)
  C = Covering(vcat(patches(C1), patches(C2)))
  for (X, Y) in keys(glueings(C1))
    add_glueing!(C, C1[X, Y])
  end
  for (X, Y) in keys(glueings(C2))
    add_glueing!(C, C2[X, Y])
  end
  return C
end

