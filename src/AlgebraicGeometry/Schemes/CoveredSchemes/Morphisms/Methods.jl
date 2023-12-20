
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
    cfg = compose(cf, cg)
    return CoveredSchemeMorphism(X, Z, cfg)
  elseif is_refinement(codomain(cf), domain(cg))[1]
    ref = refinement_morphism(codomain(cf), domain(cg))
    cf_ref_cg = compose(cf, compose(ref, cg))
    return CoveredSchemeMorphism(X, Z, cf_ref_cg, check=true)
  elseif is_refinement(domain(cg), codomain(cf))[1]
    ref = refinement_morphism(domain(cg), codomain(cf))
    A = domain(cf)
    B = codomain(cf)
    C = domain(ref)
    AxC, to_A, to_C = fiber_product(cf, ref)
    return CoveredSchemeMorphism(X, Z, compose(to_C, cg), check=true)
  else
    A = domain(cf)
    B = codomain(cf)
    C = domain(cg)
    D = codomain(cg)
    ref_B = refinement_morphism(B, default_covering(Y))
    ref_C = refinement_morphism(C, default_covering(Y))
    BxC, to_B, to_C = fiber_product(ref_B, ref_C)
    AxBxC, to_A, to_BxC = fiber_product(fc, to_B)
    return CoveredSchemeMorphism(X, Z, compose(to_BxC, compose(to_C, gc)), check=true)
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
  if get(io, :supercompact, false)
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
  if CA === CB
    AxB, to_A, to_B = fiber_product(f_cov, g_cov)
    XxY = CoveredScheme(AxB)
    to_X = CoveredSchemeMorphism(XxY, X, to_A)
    to_Y = CoveredSchemeMorphism(XxY, Y, to_B)
    return XxY, to_X, to_Y
  elseif is_refinement(CA, CB)
    inc = refinement_morphism(CA, CB)
    f_cov_inc = compose(f_cov, inc)
    AxB, to_A, to_B = fiber_product(f_cov_inc, g_cov)
    XxY = CoveredScheme(AxB)
    to_X = CoveredSchemeMorphism(XxY, X, to_A)
    to_Y = CoveredSchemeMorphism(XxY, Y, to_B)
    return XxY, to_X, to_Y
  elseif is_refinement(CB, CA)
    inc = refinement_morphism(CB, CA)
    g_cov_inc = compose(g_cov, inc)
    AxB, to_A, to_B = fiber_product(f_cov, g_cov_inc)
    XxY = CoveredScheme(AxB)
    to_X = CoveredSchemeMorphism(XxY, X, to_A)
    to_Y = CoveredSchemeMorphism(XxY, Y, to_B)
    return XxY, to_X, to_Y
  else
    C = default_covering(Z)
    inc_A = refinement_morphism(CA, C)
    inc_B = refinement_morphism(CB, C)
    CAB, to_CA, to_CB = fiber_product(inc_A, inc_B)
    AA, AA_to_A, AA_to_CAB = fiber_product(f_cov, to_CA)
    BB, BB_to_B, BB_to_CAB = fiber_product(g_cov, to_CB)
    AAxBB, to_AA, to_BB = fiber_product(AA_to_CAB, BB_to_CAB)
    XxY = CoveredScheme(AAxBB)
    to_X = CoveredSchemeMorphism(XxY, X, compose(to_AA, AA_to_A))
    to_Y = CoveredSchemeMorphism(XxY, Y, compose(to_BB, BB_to_B))
    return XxY, to_X, to_Y
  end
end

