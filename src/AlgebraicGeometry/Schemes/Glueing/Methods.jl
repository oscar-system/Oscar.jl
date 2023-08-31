########################################################################
# Printing                                                             #
########################################################################

# As for coverings and covered schemes, we need to manage offsets when printings
# coordinates so that on each open subset, the coordinates and the respective
# descriptions are all aligned on the left, in an appropriate alignment and
# distance. Moreover, since the glueing is given by some maps, we take care that
# the arrows are aligned and so arguments are aligns on the right and their
# respective images are aligned on the left.
function Base.show(io::IO, ::MIME"text/plain", G::AbsGlueing)
  io = pretty(io)
  println(io, "Glueing")
  println(io, Indent(), "of  ", Lowercase(), patches(G)[1])
  println(io, "and ", Lowercase(), patches(G)[2])
  println(io, Dedent(), "along the open subsets")
  f = glueing_morphisms(G)[1]
  pf = pullback(f)
  print(io, Indent())
  co_str = String[]
  co1 = ambient_coordinates(domain(f))
  co2 = ambient_coordinates(codomain(f))
  str = "["*join(co1, ", ")*"]"
  k1 = length(str)
  push!(co_str, str)
  str = "["*join(co2, ", ")*"]"
  k2 = length(str)
  push!(co_str, str)
  k = max(k1, k2)
  println(io, co_str[1], " "^(k-k1+3), Lowercase(), domain(f))
  print(io, co_str[2], " "^(k-k2+3), Lowercase(), codomain(f))
  if codomain(f) isa SpecOpen
    mop = maps_on_patches(f)
    if length(mop) > 0
      println(io)
      print(io, Dedent(), "defined by the map")
      length(mop) > 1 && print(io, "s")
      print(io, Indent())
      for i in 1:length(mop)
        println(io, Lowercase())
        Base.show(io, MIME"text/plain"(), mop[i])
        if i != length(mop)
          println(io)
          print(io, "------------------------------------------------------------")
        end
      end
    end
    print(io, Dedent())
  else
    c = coordinates(codomain(f))
    if length(c) > 0
      co_str = String["$(cc)" for cc in c]
      k = max(length.(co_str)...)
      println(io)
      print(io, Dedent(), "given by the pullback function")
      print(io, Indent())
      for i in 1:length(c)
        ki = length(co_str[i])
        println(io)
        print(io, " "^(k-ki)*"$(c[i]) -> $(pf(c[i]))")
      end
    end
    print(io, Dedent())
  end
end

function Base.show(io::IO, G::AbsGlueing)
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, "Glueing")
  else
    print(io, "Glueing: ", Lowercase(), patches(G)[1], " -> ", Lowercase(), patches(G)[2])
  end
end

########################################################################
# Composition of Glueings                                              #
########################################################################
@doc raw"""
    compose(G::AbsGlueing, H::AbsGlueing)

Given glueings `X ↩ U ≅ V ↪  Y` and `Y ↩ V' ≅ W ↪ Z`, return the glueing
`X ↩  V ∩ V' ↪ Z`. 

**WARNING:** In general such a glueing will not provide a separated scheme. 
Use `maximal_extension` to extend the glueing.
"""
function compose(G::AbsGlueing, H::AbsGlueing) 
  # In the default case, we rely on some conversion to the most complicated 
  # glueing and rely on the implementation for composition there.
  return compose(Glueing(G), Glueing(H))
end

function compose(G::Glueing, H::SimpleGlueing) 
  return compose(G, Glueing(H))
end

function compose(G::SimpleGlueing, H::Glueing) 
  return compose(Glueing(G), H)
end

function compose(G::Glueing, H::Glueing) 
  # make sure that Y is the second patch of the first glueing and 
  # the first patch of the second
  if patches(G)[2] === patches(H)[2]
    return compose(G, inverse(H))
  elseif patches(G)[1] === patches(H)[1]
    return compose(inverse(G), H)
  elseif patches(G)[1] === patches(H)[2]
    return compose(inverse(G), inverse(H))
  end
  X, Y = patches(G)
  Y === patches(H)[1] || error("Glueings not compatible")
  Z = patches(H)[2]
  f, f_inv = glueing_morphisms(G)
  g, g_inv = glueing_morphisms(H)
  U_new = preimage(f, domain(g))
  W_new = preimage(g_inv, codomain(f))
  V_new = intersect(codomain(f), domain(g))
  return Glueing(X, Z, 
             compose(restrict(f, U_new, V_new), restrict(g, V_new, W_new)),
             compose(restrict(g_inv, W_new, V_new), restrict(f_inv, V_new, U_new))
	     )
end

@doc raw"""
    maximal_extension(G::Glueing)

Given a glueing `X ↩ U ≅ V ↪ Y`, try to find the maximal extension to an open 
subset `U' ⊃ U` in `X` and `V' ⊃ V` in `Y` so that the resulting scheme is separated.
"""
function maximal_extension(G::Glueing)
  X = patches(G)[1]
  Y = patches(G)[2]
  f, g = glueing_morphisms(G)
  f_ext = maximal_extension(X, Y, generic_fractions(f))
  g_ext = maximal_extension(Y, X, generic_fractions(g))
  U_new = preimage(f_ext, domain(g_ext))
  V_new = preimage(g_ext, U_new)
  issubset(domain(g_ext), V_new) || error("extension failed")
  f_ext = restrict(f_ext, U_new, V_new)
  g_ext = restrict(g_ext, V_new, U_new)
  return Glueing(X, Y, f_ext, g_ext)
end

function compose(G::SimpleGlueing, H::SimpleGlueing)
  if patches(G)[2] === patches(H)[2]
    return compose(G, inverse(H))
  elseif patches(G)[1] === patches(H)[1]
    return compose(inverse(G), H)
  elseif patches(G)[1] === patches(H)[2]
    return compose(inverse(G), inverse(H))
  end
  X, Y = patches(G)
  Y === patches(H)[1] || error("Glueings not compatible")
  Z = patches(H)[2]
  f, f_inv = glueing_morphisms(G)
  g, g_inv = glueing_morphisms(H)
  U_new = PrincipalOpenSubset(ambient_scheme(domain(f)), 
                              [complement_equation(domain(f)), 
                               lifted_numerator(pullback(f)(complement_equation(domain(g))))
                              ])
  W_new = PrincipalOpenSubset(ambient_scheme(domain(g_inv)), 
                              [complement_equation(domain(g_inv)), 
                               lifted_numerator(pullback(g_inv)(complement_equation(domain(f_inv))))
                              ])
  V_new = PrincipalOpenSubset(ambient_scheme(domain(g)), 
                              [complement_equation(domain(g)), complement_equation(domain(f_inv))]
                             )
  h = compose(restrict(f, U_new, V_new, check=false), 
              restrict(g, V_new, W_new, check=false))
  h_inv = compose(restrict(g_inv, W_new, V_new, check=false),
                  restrict(f_inv, V_new, U_new, check=false))
  set_attribute!(h, :inverse, h_inv)
  set_attribute!(h_inv, :inverse, h)
  return SimpleGlueing(X, Z, h, h_inv)
end

########################################################################
# Equality tests                                                       #
########################################################################
function ==(G::AbsGlueing, H::AbsGlueing)
  if patches(G)[1] != patches(H)[1]
    return G == inverse(H)
  end
  patches(G)[2] === patches(H)[2] || return false
  glueing_morphisms(G) == glueing_morphisms(H) || return false
  return true
end

########################################################################
# Base change
########################################################################
struct BaseChangeGlueingData{T1, T2}
  phi::T1
  G::T2
  patch_change1::AbsSpecMor
  patch_change2::AbsSpecMor
end

function _compute_glueing_base_change(gd::BaseChangeGlueingData)
  # extraction of the glueing data
  G = gd.G
  pc1 = gd.patch_change1
  pc2 = gd.patch_change2
  phi = gd.phi

  # computation of the new glueing
  (X, Y) = patches(G)
  (U, V) = glueing_domains(G)
  (f, g) = glueing_morphisms(G)
  
  XX = domain(pc1)
  YY = domain(pc2)

  UU, map_U = base_change(phi, U, ambient_map=pc1)
  VV, map_V = base_change(phi, V, ambient_map=pc2)

  # TODO: The following will not yet work for glueings along SpecOpens.
  U isa PrincipalOpenSubset && V isa PrincipalOpenSubset || error("base change not implemented for glueings along SpecOpens")

  _, ff, _ = base_change(phi, f, domain_map=map_U, codomain_map=map_V)
  _, gg, _ = base_change(phi, g, domain_map=map_V, codomain_map=map_U)

  return Glueing(XX, YY, ff, gg)
end

function base_change(phi::Any, G::AbsGlueing;
    patch_change1::AbsSpecMor=base_change(phi, glueing_domains(G)[1])[2], # The base change morphism for 
    patch_change2::AbsSpecMor=base_change(phi, glueing_domains(G)[2])[2]  # the two glueing patches
  )
  gd = BaseChangeGlueingData{typeof(phi), typeof(G)}(phi, G, patch_change1, patch_change2)
  return LazyGlueing(domain(patch_change1), domain(patch_change2), _compute_glueing_base_change, gd)
end

