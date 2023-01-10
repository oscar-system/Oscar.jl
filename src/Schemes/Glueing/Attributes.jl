export underlying_glueing, patches, glueing_morphisms, inverse, glueing_domains

########################################################################
# The interface for abstract glueings of affine schemes                #
########################################################################
function underlying_glueing(G::AbsGlueing)
  error("trying to call for `underlying_glueing` of $G but nothing is implemented")
end

@Markdown.doc """
    patches(G::AbsGlueing)

Return a pair of affine schemes `(X, Y)` which are glued by `G`.
"""
function patches(G::AbsGlueing)
  return patches(underlying_glueing(G))
end

@Markdown.doc """
    glueing_morphisms(G::AbsGlueing)

Return a pair of mutually inverse isomorphisms `(f, g)` 
of open subsets ``U`` and ``V`` of the respective 
`patches` of `G` which are used for the glueing 
identification.
"""
function glueing_morphisms(G::AbsGlueing)
  return glueing_morphisms(underlying_glueing(G))
end

@Markdown.doc """
    glueing_domains(G::AbsGlueing)

Return a pair of open subsets ``U`` and ``V`` of the 
respective `patches` of `G` which are glued by `G` 
along the `glueing_morphisms` of `G`.
"""
function glueing_domains(G::AbsGlueing)
  return glueing_domains(underlying_glueing(G))
end

@Markdown.doc """
    inverse(G::AbsGlueing)

Return the glueing `H` with `patches`, `glueing_domains`, 
and `glueing_morphisms` in opposite order compared to `G`.
"""
@attr function inverse(G::AbsGlueing)
  Ginv = inverse(underlying_glueing(G))
  set_attribute!(Ginv, :inverse, G)
  set_attribute!(G, :inverse, Ginv)
  return Ginv
end

########################################################################
# Getters for Glueing                                                  #
########################################################################
patches(G::Glueing) = G.X, G.Y
glueing_morphisms(G::Glueing) = G.f, G.g
glueing_domains(G::Glueing) = domain(G.f), domain(G.g)
@attr function inverse(G::Glueing)
  Ginv = Glueing(G.Y, G.X, G.g, G.f, check=false)
  set_attribute!(Ginv, :inverse, G)
  return Ginv
end

########################################################################
# Getters for SimpleGlueing                                            #
########################################################################
patches(G::SimpleGlueing) = (G.X, G.Y)
glueing_morphisms(G::SimpleGlueing) = (G.f, G.g)
glueing_domains(G::SimpleGlueing) = (G.U, G.V)

@attr SimpleGlueing function inverse(G::SimpleGlueing)
  Ginv = SimpleGlueing(G.Y, G.X, G.g, G.f, check=false)
  set_attribute!(Ginv, :inverse, G)
  return Ginv
end


