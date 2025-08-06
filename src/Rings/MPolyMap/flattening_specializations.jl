# specialized routines with improvements

# some unifying getters to allow for generic code
_poly_lift(x::MPolyRingElem) = x
_poly_lift(x::MPolyQuoRingElem) = x.f

_poly_ring(R::MPolyRing) = R
_poly_ring(R::MPolyQuoRing) = base_ring(R)

function (phi::RingFlattening{<:Union{
                                    <:MPolyRing{<:MPolyRingElem},
                                    <:MPolyRing{<:MPolyQuoRingElem},
                                    <:MPolyQuoRing{<:MPolyRingElem{<:MPolyRingElem}},
                                    <:MPolyQuoRing{<:MPolyRingElem{<:MPolyQuoRingElem}}
                                   }, 
                              <:Union{
                                    <:MPolyRing, <:MPolyQuoRing
                                   }})(x::RingElem)
  !phi.cached && return _compute_flattened_element(phi, x)
  return get!(flat_counterparts(phi), x) do
    _compute_flattened_element(phi, x)
  end::elem_type(codomain(phi))
end

function _compute_flattened_element(phi, x)
  is_zero(x) && return zero(codomain(phi))
  R = domain(phi)
  parent(x) === R || return _compute_flattened_element(phi, R(x))
  RP = _poly_ring(R)
  P = coefficient_ring(RP)
  PP = _poly_ring(P)
  kk = coefficient_ring(PP)
  S = codomain(phi)
  SP = _poly_ring(S)
  @assert kk === coefficient_ring(SP)
  if is_constant(_poly_lift(x))
    c = first(AbstractAlgebra.coefficients(_poly_lift(x)))
    is_constant(_poly_lift(c)) && return S(first(AbstractAlgebra.coefficients(_poly_lift(c))))
  end
  res_ctx = MPolyBuildCtx(SP)
  for (c, e) in zip(AbstractAlgebra.coefficients(_poly_lift(x)), AbstractAlgebra.exponent_vectors(_poly_lift(x)))
    for (cc, ee) in zip(AbstractAlgebra.coefficients(_poly_lift(c)), AbstractAlgebra.exponent_vectors(_poly_lift(c)))
      push_term!(res_ctx, cc, vcat(e, ee))
    end
  end
  return S(finish(res_ctx))
end

# TODO: One can still try to do something like the above optimizations 
# also for polynomial rings over localizations of polynomial rings. The 
# problem there is that one can not directly concatenate exponent vectors 
# because we have potential denominators floating around. But should this 
# turn out to be a bottleneck, one can probably think of something.

# TODO: Think of a way to speed up the inverse of a RingFlattening.
# The problem is that at the moment the inverse is a standalone map, 
# stored in the internals of a RingFlattening. When it's applied to an 
# element, it does not know about its origin anymore, so we can not 
# overwrite this via some dispatch for RingFlattening. 
