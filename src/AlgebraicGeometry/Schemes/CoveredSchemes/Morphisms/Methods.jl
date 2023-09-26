
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
    print(io, "Morphism")
  else
    print(io, "Morphism: ", Lowercase(), domain(f), " -> ", Lowercase(), codomain(f))
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
  println(io, "Morphism")
  print(io, Indent(), "from ", Lowercase())
  Oscar._show_semi_compact(io, domain(f), domain(g), "a")
  println(io)
  print(io, "to   ", Lowercase())
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

function fiber_product(
    i1::CoveredClosedEmbedding,
    i2::CoveredClosedEmbedding
  )
  X1 = domain(i1)
  X2 = domain(i2)
  Y = codomain(i1)
  Y === codomain(i2) || error("codomains do not coincide")
  i1_cov = covering_morphism(i1)
  i2_cov = covering_morphism(i2)
  codomain(i1_cov) === codomain(i2_cov) || error("case of different coverings in codomain not implemented")
  cod_cov = codomain(i1_cov)
  #=
  cod_ref, ref1, ref2 = common_refinement(codomain(i1_cov), codomain(i2_cov))
  dom_ref1, i1_res = fiber_product(i1, ref1)
  dom_ref2, i2_res = fiber_product(i1, ref2)
  # etc. etc.... This is roughly the generic code to come.
  =#
  I1 = image_ideal(i1)
  pb_I1 = pullback(i2, I1)
  I2 = image_ideal(i2)
  pb_I2 = pullback(i1, I2)

  j1 = Oscar.CoveredClosedEmbedding(domain(i2), pb_I1)
  Z = domain(j1)
  morphism_dict = IdDict{AbsSpec, ClosedEmbedding}()
  for U in affine_charts(Z)
    V2 = codomain(j1[U])
    W = codomain(i2[V2])
    V1_candidates = maps_with_given_codomain(i1, W)
    @assert length(V1_candidates) == 1 "not the correct number of patches found"
    V1 = domain(first(V1_candidates))
    x = gens(OO(V1))
    @show x
    lift_x = [preimage(pullback(i1[V1]), f) for f in x]
    @show lift_x
    pb_x = pullback(i2[V2]).(lift_x)
    @show pb_x
    pb_x = pullback(j1[U]).(pb_x)
    @show pb_x
    morphism_dict[U] = ClosedEmbedding(SpecMor(U, V1, pb_x, check=false), pb_I2(V1), check=false)
  end
  j2_cov = CoveringMorphism(default_covering(Z), domain(i1_cov), morphism_dict, check=false)
  j2 = Oscar.CoveredClosedEmbedding(Z, X1, j2_cov)
  return j1, j2
end
