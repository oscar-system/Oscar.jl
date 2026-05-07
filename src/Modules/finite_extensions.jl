underlying_map(psi::FiniteExtension) = psi.phi
domain(psi::FiniteExtension{DT}) where DT = domain(underlying_map(psi))::DT
codomain(psi::FiniteExtension{DT, CT}) where {DT, CT} = codomain(underlying_map(psi))::CT

function fiber_ideal(psi::FiniteExtension)
  if !isdefined(psi, :fiber_ideal)
    R = codomain(psi)
    psi.fiber_ideal = ideal(R, images_of_generators(psi.phi))
  end
  return psi.fiber_ideal::ideal_type(codomain(psi))
end

function fiber_quotient(psi::FiniteExtension)
  if !isdefined(psi, :fiber_quotient)
    Q, _ = quo(codomain(psi), fiber_ideal(psi))
    psi.fiber_quotient = Q
  end
  return psi.fiber_quotient
end

function basis(psi::FiniteExtension)
  if !isdefined(psi, :basis)
    codomain_as_module(psi) # fill the cache
    # R = codomain(psi)
    # kk = coefficient_ring(R)
    # psi.basis = vector_space_basis(kk, fiber_quotient(psi))
  end
  return psi.basis::Vector{elem_type(base_ring(codomain(psi)))}
end

function vector_space_basis(kk::Field, Q::MPolyQuoRing)
  @assert kk === coefficient_ring(Q)
  I = modulus(Q)
  P = base_ring(I)
  SI = singular_generators(groebner_basis(I))
  sb = Singular.kbase(SI)
  return [P(x) for x in gens(sb)]::Vector{elem_type(P)}
end

function codomain_as_module(psi::FiniteExtension)
  if !isdefined(psi, :codomain_as_module)
    phi = underlying_map(psi)
    g, rel, interp = present_finite_extension_ring(phi)
    B = domain(psi)
    R = codomain(psi)
    P = base_ring(R)
    G = grading_group(B)
    F = graded_free_module(B, [degree(Int, x) for x in g])
    rels = isempty(rel) ? elem_type(F)[] : elem_type(F)[sum(rel[i, j]*F[j] for j in 1:ncols(rel); init=zero(F)) for i in 1:nrows(rel)]
    psi.codomain_as_module = SubquoModule(F, gens(F), rels)
    psi.interp = interp
    psi.basis = g
  end
  return psi.codomain_as_module::SubquoModule{elem_type(domain(psi))}
end

function pushforward(psi::FiniteExtension, F::FreeMod)
  M = codomain_as_module(psi)
  summands = [twist(M, Int(g[1])) for g in degrees_of_generators(F)]
  result, _ = direct_sum(summands...)
  interp = function(v::FreeModElem)
    res = zero(result)
    for (i, c) in coordinates(v)
      inj = canonical_injection(result, i)
      w = psi.interp(lift(c))
      ww = sum(cc*M[j] for (j, cc) in enumerate(w); init=zero(M))
      res += inj(summands[i](coordinates(ww)))
    end
    return res
  end
  interp_inv = function(w::SubquoModuleElem)
    g = basis(psi)
    R = codomain(psi)
    res = zero(F)
    phi = underlying_map(psi)
    for i in 1:length(summands)
      pr = canonical_projection(result, i)
      ww = pr(w)
      cw = coordinates(ww)
      res = sum(phi(c)*g[j]*F[i] for (j, c) in cw; init=res)
    end
    return res
  end
  return result, MapFromFunc(F, result, interp, interp_inv)
end

function pushforward(psi::FiniteExtension, M::SubquoModule;
    ambient_pushforward=pushforward(psi, ambient_free_module(M))
  )
  psiF, interp = ambient_pushforward
  psi_gens = elem_type(psiF)[interp(g*v) for v in ambient_representatives_generators(M) for g in basis(psi)]
  psi_rels = elem_type(psiF)[interp(g*v) for v in relations(M) for g in basis(psi)]
  I, _ = sub(psiF, psi_gens)
  res, _ = quo(I, psi_rels)
  return res
end

