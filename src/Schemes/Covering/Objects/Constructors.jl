export empty_covering, disjoint_union, standard_covering
########################################################################
# Constructors for Covering                                            #
########################################################################
### The default constructor
# Returns a scheme in which every affine patch is only 
# glued to itself via the identity.
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

