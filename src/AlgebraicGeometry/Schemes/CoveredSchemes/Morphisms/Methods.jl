
########################################################################
# Methods for CoveredSchemeMorphism                                    #
########################################################################
@doc raw"""
    simplify!(X::AbsCoveredScheme)

Apply `simplify` to the `default_covering` of `X` and store the 
resulting `Covering` ``C'`` and its identification with the default 
covering in `X`; return a tuple ``(X, C')``.
"""
function simplify!(X::AbsCoveredScheme)
  C = default_covering(X)
  Csimp, i, j = simplify(C)
  push!(coverings(X), Csimp)
  refinements(X)[(C, Csimp)] = j
  refinements(X)[(Csimp, C)] = i
  set_attribute!(X, :simplified_covering, Csimp)
  return X
end


########################################################################
# Auxiliary methods for compatibility                                  #
########################################################################
lifted_numerator(f::MPolyRingElem) = f
lifted_denominator(f::MPolyRingElem) = one(f)
lifted_numerator(f::MPolyQuoRingElem) = lift(f)
lifted_denominator(f::MPolyQuoRingElem) = one(lift(f))

function compose(f::AbsCoveredSchemeMorphism, g::AbsCoveredSchemeMorphism)
  X = domain(f)
  Y = codomain(f)
  Y === domain(g) || error("maps not compatible")
  Z = codomain(g)
  cf = covering_morphism(f)
  cg = covering_morphism(g)
  if codomain(cf) === domain(cg) 
    mor_dict = IdDict{AbsSpec, AbsSpecMor}() # TODO: Keep the type of the morphisms?
    for U in basic_patches(domain(cf))
      mor_dict[U] = compose(cf[U], cg[codomain(cf[U])])
    end
    cc = CoveringMorphism(domain(cf), codomain(cg), mor_dict, check=false)
    return CoveredSchemeMorphism(X, Z, cc)
  else
    haskey(refinements(Y), (codomain(cf), domain(cg))) || error("composition of this complicated case is not yet implemented")
    ref_mor = Y[codomain(cf), domain(cg)] # the refinement `CoveringMorphism`
    mor_dict = IdDict{AbsSpec, AbsSpecMor}() # TODO: Keep the type of the morphisms?
    for U in basic_patches(domain(cf))
      V = codomain(cf[U])
      W = codomain(ref_mor[V])
      mor_dict[U] = compose(compose(cf[U], ref_mor[V]), cg[W])
    end
    cc = CoveringMorphism(domain(cf), codomain(cg), mor_dict, check=false)
    return CoveredSchemeMorphism(X, Z, cc)
  end
end

function maps_with_given_codomain(f::AbsCoveredSchemeMorphism, V::AbsSpec)
  fcov = covering_morphism(f)
  result = Vector{AbsSpecMor}()
  for U in keys(morphisms(fcov))
    floc = morphisms(fcov)[U]
    codomain(floc) === V || continue
    push!(result, floc)
  end
  return result
end

########################################################################
# Comparison
########################################################################
function ==(f::AbsCoveredSchemeMorphism, g::AbsCoveredSchemeMorphism)
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  f_cov = covering_morphism(f)
  g_cov = covering_morphism(g)
  domain(f_cov) === domain(g_cov) || error("comparison across refinements not implemented")
  codomain(f_cov) === codomain(g_cov) || error("comparison across refinements not implemented")
  return all(U->(f_cov[U] == g_cov[U]), patches(domain(f_cov)))
end

########################################################################
# Base change
########################################################################
function base_change(phi::Any, f::AbsCoveredSchemeMorphism;
    domain_map::AbsCoveredSchemeMorphism=base_change(phi, domain(f))[2],
    codomain_map::AbsCoveredSchemeMorphism=base_change(phi, codomain(f))[2]
  )
  f_cov = covering_morphism(f)
  dom_cov = covering_morphism(domain_map)
  cod_cov = covering_morphism(codomain_map)
  _, ff_cov_map, _ = base_change(phi, f_cov, domain_map=dom_cov, codomain_map=cod_cov)
  X = domain(f)
  Y = codomain(f)
  XX = domain(domain_map)
  YY = domain(codomain_map)
  return domain_map, CoveredSchemeMorphism(XX, YY, ff_cov_map), codomain_map
end

function _register_birationality!(f::AbsCoveredSchemeMorphism, 
    g::AbsSpecMor, ginv::AbsSpecMor)
  set_attribute!(g, :inverse, ginv)
  set_attribute!(ginv, :inverse, g)
  return _register_birationality(f, g)
end

function _register_birationality!(f::AbsCoveredSchemeMorphism, 
    g::AbsSpecMor
  )
  set_attribute!(f, :is_birational, true)
  set_attribute!(f, :iso_on_open_subset, g)
end

###############################################################################
#
#  Printing
#
###############################################################################

# We use a pattern for printings morphisms, glueings, etc...
#
# In supercompact printing, we just write what it is, super shortly.
# For normal compact printing, we mention what it is, then use colons to
# describe "domain -> codomain".
function Base.show(io::IO, f::AbsCoveredSchemeMorphism)
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, "Covered scheme morphism")
  else
    print(io, "Hom: ", Lowercase(), domain(f), " -> ", Lowercase(), codomain(f))
  end
end

# Here the `_show_semi_compact` allows us to avoid the redundancy on the
# printing of the domain/codomain, choose the covering of the scheme to print,
# and associate a letter to the labels of the charts - "a" for the domain and
# "b" for the codomain.
#
# We also have a `_show_semi_compact` for the associate covering morphism, where
# again we just avoid to re-write what are the domain and codomain.
function Base.show(io::IO, ::MIME"text/plain", f::AbsCoveredSchemeMorphism)
  io = pretty(io)
  g = covering_morphism(f)
  println(io, "Covered scheme morphism")
  print(io, Indent(), "from ", Lowercase())
  Oscar._show_semi_compact(io, domain(f), domain(g), "a")
  println(io)
  print(io, "to ", Lowercase())
  Oscar._show_semi_compact(io, codomain(f), codomain(g), "b")
  if min(length(domain(g)), length(codomain(g))) == 0
    print(io, Dedent())
  else
    println(io, Dedent())
    print(io, "given by the pullback function")
    length(domain(g)) > 1 && print(io, "s")
    println(io, Indent())
    Oscar._show_semi_compact(io, covering_morphism(f))
  end
end

