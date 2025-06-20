
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
  # The general problem here is that the `CoveringMorphism`s for `f` and 
  # `g` need not have compatible domains and codomains. If they do, we're 
  # lucky, but if they don't we first need to create the necessary 
  # refinements and then compose the induced maps.
  if codomain(cf) === domain(cg) 
    # The easy case
    cfg = compose(cf, cg)
    return CoveredSchemeMorphism(X, Z, cfg; check=false)
  elseif is_refinement(codomain(cf), domain(cg))[1]
    # Another rather easy case: We only have to put the refinement morphism 
    # in the middle.
    ref = refinement_morphism(codomain(cf), domain(cg))
    cf_ref_cg = compose(cf, compose(ref, cg))
    return CoveredSchemeMorphism(X, Z, cf_ref_cg, check=false)
  elseif is_refinement(domain(cg), codomain(cf))[1]
    # A more tricky case: 
    #
    # AxB' --------> B' ---cg---> C
    #  |             |
    #  V             V ref
    #  A --cf------> B
    #
    # Here we first have to complete the square to get another 
    # `CoveringMorphism` representing `f` with `B'` as is codomain.
    ref = refinement_morphism(domain(cg), codomain(cf))
    A = domain(cf)
    B = codomain(cf)
    C = domain(ref)
    AxC, to_A, to_C = fiber_product(cf, ref)
    return CoveredSchemeMorphism(X, Z, compose(to_C, cg), check=false)
  else
    # The most complicated case:
    #
    #  AxBxC ----> BxC
    #   /          / \
    #  V          V   V
    # A ---cf--> B     C ---cg--> D
    #             \   /
    #              V V
    #         default_covering(Y)
    #
    # Here we have to complete two squares, first (B->def, C->def)
    # and then (A->B, BxC->B).
    A = domain(cf)
    B = codomain(cf)
    C = domain(cg)
    D = codomain(cg)
    ref_B = refinement_morphism(B, default_covering(Y))
    ref_C = refinement_morphism(C, default_covering(Y))
    BxC, to_B, to_C = fiber_product(ref_B, ref_C)
    AxBxC, to_A, to_BxC = fiber_product(cf, to_B)
    return CoveredSchemeMorphism(X, Z, compose(to_BxC, compose(to_C, cg)), check=false)
  end
end

function maps_with_given_codomain(f::AbsCoveredSchemeMorphism, V::AbsAffineScheme)
  fcov = covering_morphism(f)
  result = Vector{AbsAffineSchemeMor}()
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
  @assert codomain(domain_map) === domain(f) "domains do not match"
  @assert codomain(codomain_map) === codomain(f) "codomains do not match"
  f_cov = covering_morphism(f)
  dom_cov = covering_morphism(domain_map)
  cod_cov = covering_morphism(codomain_map)
  # The codomains of the base change maps need not be compatible with 
  # the domain/codomain of f_cov. If this is the case, we need to refine 
  # all these.
  if codomain(dom_cov) !== domain(f_cov)
    if is_refinement(codomain(dom_cov), domain(f_cov))[1]
      ref = refinement_morphism(codomain(dom_cov), domain(f_cov))
      dom_cov = compose(dom_cov, ref)
    elseif is_refinement(domain(f_cov), codomain(dom_cov))
      # This should be the most common case, really: base_change 
      # has happened on the `default_covering`s and `f` has a refinement 
      # thereof as its codomain.
      ref = refinement_morphism(domain(f_cov), codomain(dom_cov))
      _, to_dom_dom_cov, to_dom_f_cov = fiber_product(dom_cov, ref)
      dom_cov = to_dom_f_cov
    else
      # This should not really happen usually, since we assume a base_change 
      # to be carried out on the default_covering.
      error("case not implemented")
    end
  end
  if codomain(cod_cov) !== codomain(f_cov)
    # We must assume that `cod_cov` is realized w.r.t. the `default_covering`s 
    # on both sides. Otherwise, we have no chance to write down the lifting 
    # map to the rings with the new coefficient ring.
    ref = refinement_morphism(codomain(f_cov), default_covering(codomain(f)))
    f_cov = compose(f_cov, ref)
  end
    
  @assert codomain(dom_cov) === domain(f_cov)
  @assert codomain(cod_cov) === codomain(f_cov) "base change in the codomain is not possible unless one is using the `default_covering`"
  _, ff_cov_map, _ = base_change(phi, f_cov, domain_map=dom_cov, codomain_map=cod_cov)
  X = domain(f)
  Y = codomain(f)
  XX = domain(domain_map)
  YY = domain(codomain_map)
  return domain_map, CoveredSchemeMorphism(XX, YY, ff_cov_map; check=false), codomain_map
end

function _register_birationality!(f::AbsCoveredSchemeMorphism, 
    g::AbsAffineSchemeMor, ginv::AbsAffineSchemeMor)
  set_attribute!(g, :inverse, ginv)
  set_attribute!(ginv, :inverse, g)
  return _register_birationality(f, g)
end

function _register_birationality!(f::AbsCoveredSchemeMorphism, 
    g::AbsAffineSchemeMor
  )
  set_attribute!(f, :is_birational, true)
  set_attribute!(f, :iso_on_open_subset, g)
end

###############################################################################
#
#  Printing
#
###############################################################################

# We use a pattern for printings morphisms, gluings, etc...
#
# In terse printing, we just write what it is, super shortly.
# For normal compact printing, we mention what it is, then use colons to
# describe "domain -> codomain".
function Base.show(io::IO, f::AbsCoveredSchemeMorphism)
  if is_terse(io)
    print(io, "Covered scheme morphism")
  else
    io = pretty(io)
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
  show(IOContext(io, :show_semi_compact => true, :covering => domain(g), :label => "a"), domain(f))
  println(io)
  print(io, "to ", Lowercase())
  show(IOContext(io, :show_semi_compact => true, :covering => codomain(g), :label => "b"), codomain(f))
  if min(length(domain(g)), length(codomain(g))) == 0
    print(io, Dedent())
  else
    println(io, Dedent())
    print(io, "given by the pullback function")
    length(domain(g)) != 1 && print(io, "s")
    println(io, Indent())
    show(IOContext(io, :show_semi_compact => true), covering_morphism(f))
  end
end

########################################################################
# fiber products
########################################################################

@doc raw"""
    fiber_product(f::AbsCoveredSchemeMorphism, g::AbsCoveredSchemeMorphism)

For a diagram 
    XxY ----> Y
     |        | g
     V        V
     X------> Z
         f
this computes the fiber product `XxY` together with the canonical maps 
to `X` and `Y` and returns the resulting triple.
"""
function fiber_product(f::AbsCoveredSchemeMorphism, g::AbsCoveredSchemeMorphism)
  X = domain(f)
  Y = domain(g)
  Z = codomain(f)
  @assert Z === codomain(g)
  f_cov = covering_morphism(f)
  g_cov = covering_morphism(g)
  A = domain(f_cov)
  B = domain(g_cov)
  CA = codomain(f_cov)
  CB = codomain(g_cov)
  # The problem is that the `CoveringMorphism`s representing `f` and `g` need 
  # not have compatible codomains. Hence we will need to pass to necessary 
  # refinements and their induced maps in most cases.
  if CA === CB
    # The easy case.
    AxB, to_A, to_B = fiber_product(f_cov, g_cov)
    XxY = CoveredScheme(AxB)
    to_X = CoveredSchemeMorphism(XxY, X, to_A; check=false)
    to_Y = CoveredSchemeMorphism(XxY, Y, to_B; check=false)
    return XxY, to_X, to_Y
  elseif is_refinement(CA, CB)[1]
    # We have to complete the square
    #
    #                              B
    #                              |
    #                             g_cov
    #                              |
    #                              V
    # A ---f_cov---> CA ---ref---> CB
    #
    # which is still rather easy.
    inc = refinement_morphism(CA, CB)
    f_cov_inc = compose(f_cov, inc)
    AxB, to_A, to_B = fiber_product(f_cov_inc, g_cov)
    XxY = CoveredScheme(AxB)
    to_X = CoveredSchemeMorphism(XxY, X, to_A; check=false)
    to_Y = CoveredSchemeMorphism(XxY, Y, to_B; check=false)
    return XxY, to_X, to_Y
  elseif is_refinement(CB, CA)[1]
    # Similar to the above case
    inc = refinement_morphism(CB, CA)
    g_cov_inc = compose(g_cov, inc)
    AxB, to_A, to_B = fiber_product(f_cov, g_cov_inc)
    XxY = CoveredScheme(AxB)
    to_X = CoveredSchemeMorphism(XxY, X, to_A; check=false)
    to_Y = CoveredSchemeMorphism(XxY, Y, to_B; check=false)
    return XxY, to_X, to_Y
  else
    # In this case we complete the following square
    # successively:
    #
    # AxCAxCBxB---> CAxCBxB------> B
    #  |             |             |
    #  |             |            g_cov
    #  |             |             |
    #  V             V             V
    # AxCAxCB-----> CAxCB -------> CB
    # |              |             |
    # |              |            inc_B
    # |              |             |
    # V              V             V
    # A ---f_cov---> CA ---inc_A-> C 
    #
    # TODO: Maybe it's easier to first compose the maps 
    # on the boundary and do only one square? We should think 
    # about it!
    C = default_covering(Z)
    inc_A = refinement_morphism(CA, C)
    inc_B = refinement_morphism(CB, C)
    CAB, to_CA, to_CB = fiber_product(inc_A, inc_B)
    AA, AA_to_A, AA_to_CAB = fiber_product(f_cov, to_CA)
    BB, BB_to_B, BB_to_CAB = fiber_product(g_cov, to_CB)
    AAxBB, to_AA, to_BB = fiber_product(AA_to_CAB, BB_to_CAB)
    XxY = CoveredScheme(AAxBB)
    to_X = CoveredSchemeMorphism(XxY, X, compose(to_AA, AA_to_A); check=false)
    to_Y = CoveredSchemeMorphism(XxY, Y, compose(to_BB, BB_to_B); check=false)
    return XxY, to_X, to_Y
  end
end
