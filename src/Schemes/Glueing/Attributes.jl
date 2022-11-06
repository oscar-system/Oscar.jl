export underlying_glueing, patches, glueing_morphisms, inverse, glueing_domains

########################################################################
# The interface for abstract glueings of affine schemes                #
########################################################################
function underlying_glueing(G::AbsGlueing)
  error("trying to call for `underlying_glueing` of $G but nothing is implemented")
end

function patches(G::AbsGlueing)
  return patches(underlying_glueing(G))
end

function glueing_morphisms(G::AbsGlueing)
  return glueing_morphisms(underlying_glueing(G))
end

function glueing_domains(G::AbsGlueing)
  return glueing_domains(underlying_glueing(G))
end

@attr function inverse(G::AbsGlueing)
  Ginv = G(inverse(underlying_glueing(G)))
  set_attribute!(Ginv, :inverse, G)
  return Ginv
end

########################################################################
# Getters for Glueing                                                  #
########################################################################
patches(G::Glueing) = G.X, G.Y
glueing_morphisms(G::Glueing) = G.f, G.g
glueing_domains(G::Glueing) = domain(G.f), domain(G.g)
inverse(G::Glueing) = Glueing(G.Y, G.X, G.g, G.f, check=false)

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


