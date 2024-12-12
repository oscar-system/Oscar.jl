# specialized routines with improvements
function (phi::RingFlattening{<:MPolyRing{<:MPolyRingElem}, <:MPolyRing})(x::MPolyRingElem)
  R = domain(phi)
  parent(x) === R || return phi(R(x))
  P = coefficient_ring(R)
  kk = coefficient_ring(P)
  S = codomain(phi)
  @assert kk === coefficient_ring(S)
  res_ctx = MPolyBuildCtx(S)
  for (c, e) in zip(AbstractAlgebra.coefficients(x), AbstractAlgebra.exponent_vectors(x))
    for (cc, ee) in zip(AbstractAlgebra.coefficients(c), AbstractAlgebra.exponent_vectors(c))
      push_term!(res_ctx, cc, vcat(e, ee))
    end
  end
  return finish(res_ctx)
end

