
########################################################################
# The interface for abstract gluings of affine schemes                #
########################################################################
function underlying_gluing(G::AbsGluing)
  error("trying to call for `underlying_gluing` of $G but nothing is implemented")
end

@doc raw"""
    patches(G::AbsGluing)

Return a pair of affine schemes `(X, Y)` which are glued by `G`.
"""
function patches(G::AbsGluing)
  return patches(underlying_gluing(G))
end

@doc raw"""
    gluing_morphisms(G::AbsGluing)

Return a pair of mutually inverse isomorphisms `(f, g)` 
of open subsets ``U`` and ``V`` of the respective 
`patches` of `G` which are used for the gluing 
identification.
"""
function gluing_morphisms(G::AbsGluing)
  return gluing_morphisms(underlying_gluing(G))
end

@doc raw"""
    gluing_domains(G::AbsGluing)

Return a pair of open subsets ``U`` and ``V`` of the 
respective `patches` of `G` which are glued by `G` 
along the `gluing_morphisms` of `G`.
"""
function gluing_domains(G::AbsGluing)
  return gluing_domains(underlying_gluing(G))
end

@doc raw"""
    inverse(G::AbsGluing)

Return the gluing `H` with `patches`, `gluing_domains`, 
and `gluing_morphisms` in opposite order compared to `G`.
"""
@attr Any function inverse(G::AbsGluing)
  Ginv = inverse(underlying_gluing(G))
  set_attribute!(Ginv, :inverse, G)
  set_attribute!(G, :inverse, Ginv)
  return Ginv
end

########################################################################
# Getters for Gluing                                                  #
########################################################################
patches(G::Gluing) = G.X, G.Y
gluing_morphisms(G::Gluing) = G.f, G.g
gluing_domains(G::Gluing) = domain(G.f), domain(G.g)
@attr Any function inverse(G::Gluing)
  Ginv = Gluing(G.Y, G.X, G.g, G.f, check=false)
  set_attribute!(Ginv, :inverse, G)
  return Ginv
end

########################################################################
# Getters for SimpleGluing                                            #
########################################################################
patches(G::SimpleGluing) = (G.X, G.Y)
gluing_morphisms(G::SimpleGluing) = (G.f, G.g)
gluing_domains(G::SimpleGluing) = (G.U, G.V)

@attr AbsGluing function inverse(G::SimpleGluing)
  Ginv = SimpleGluing(G.Y, G.X, G.g, G.f, check=false)
  set_attribute!(Ginv, :inverse, G)
  return Ginv
end


