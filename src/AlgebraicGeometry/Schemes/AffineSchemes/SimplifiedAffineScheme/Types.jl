########################################################################
# Simplified Spectra                                                   #
########################################################################
@attributes mutable struct SimplifiedAffineScheme{BaseRingType, RingType<:Ring} <: AbsAffineScheme{BaseRingType, RingType}
  X::AbsAffineScheme
  Y::AbsAffineScheme
  f::AbsAffineSchemeMor
  g::AbsAffineSchemeMor

  function SimplifiedAffineScheme(X::AbsAffineScheme, Y::AbsAffineScheme, f::AbsAffineSchemeMor, g::AbsAffineSchemeMor;
      check::Bool=true
    )
    domain(f) === X || error("map is not compatible")
    codomain(f) === Y || error("map is not compatible")
    domain(g) === Y || error("map is not compatible")
    codomain(g) === X || error("map is not compatible")

    @check is_identity_map(compose(f, g)) && is_identity_map(compose(g, f)) "maps are not inverse to each other"

    result = new{typeof(base_ring(X)), typeof(OO(X))}(X, Y)
    # We need to rewrap the identification maps so that the (co-)domains match
    fwrap = morphism(result, Y, pullback(f), check=false)
    gwrap = morphism(Y, result, pullback(g), check=false)
    set_attribute!(fwrap, :inverse, gwrap)
    set_attribute!(gwrap, :inverse, fwrap)
    result.f = fwrap
    result.g = gwrap
    return result 
  end
end
