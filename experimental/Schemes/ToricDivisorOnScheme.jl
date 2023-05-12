########################################################################
# Toric divisors on (toric) schemes
#
# This is a wrapper for the `ToricDivisor` functionality to port them 
# to the realm of `ToricCoveredScheme`s. All the functionality is lazy 
# and should by default only trigger computations on the toric side. 
# However, despite inheriting all those methods, the actual type fits 
# within the scheme's framework and the truely geometric objects are 
# derived and cached via the respective calls. In particular, these 
# divisors can seamlessly be used with usual divisors. 
########################################################################
mutable struct ToricDivisorOnScheme{CoveredSchemeType<:ToricCoveredScheme,
                                    CoeffRingType<:AbstractAlgebra.Ring, 
                                    CoeffRingElemType<:AbstractAlgebra.RingElem
                                   }<:AbsAlgebraicCycle{CoveredSchemeType, CoeffRingType}

  td::ToricDivisor # The original toric divisor
  X::CoveredSchemeType # The scheme version of the toric variety; 
                       # This must be set since otherwise arithmetic 
                       # of divisors will have massive incompatibilities!
  # Next comes the Weil divisor version of `td`. This is only computed on demand. 
  d::WeilDivisor{CoveredSchemeType, CoeffRingType, CoeffRingElemType}

  function ToricDivisorOnScheme(X::ToricCoveredScheme, td::ToricDivisor)
    return new{typeof(X), typeof(ZZ), elem_type(ZZ)}(td, X)
  end

  function ToricDivisorOnScheme(td::ToricDivisor)
    ntv = toric_variety(td)
    X = ToricCoveredScheme(ntv)
    return new{typeof(X), typeof(ZZ), elem_type(ZZ)}(td, X)
  end
end

### basic getters 
function underlying_toric_divisor(tdos::ToricDivisorOnScheme)
  return tdos.td
end

function underlying_weil_divisor(tdos::ToricDivisorOnScheme)
  if !isdefined(tdos, :d)
    # Compute the actual Weil divisor and store it in the field .d
    error("computation of associated Weil divisor not implemented")
  end
  return tdos.d
end

# forwarding the functionality for AbsAlgebraicCycle
function underlying_cycle(tdos::ToricDivisorOnScheme)
  if !isdefined(tdos, :d)
    # Compute the actual Weil divisor and store it in the field .d
    error("computation of associated Weil divisor not implemented")
  end
  return tdos.d
end

function base_scheme(tdos::ToricDivisorOnScheme)
  return tdos.X
end

########################################################################
# Inherited functionality from the toric divisors 
#
# This list was made by going through the output 
# of `methodswith(ToricDivisor)`.
########################################################################

# Because of ambiguity, we need to spell the following methods out completely.
function *(c::Integer, tdos::ToricDivisorOnScheme)
  result = ToricDivisorOnScheme(tdos.X, c*underlying_toric_divisor(tdos))
  if isdefined(tdos, :d)
    result.d = c*underlying_weil_divisor(tdos)
  end
  return result
end

function *(c::Int, tdos::ToricDivisorOnScheme)
  result = ToricDivisorOnScheme(tdos.X, c*underlying_toric_divisor(tdos))
  if isdefined(tdos, :d)
    result.d = c*underlying_weil_divisor(tdos)
  end
  return result
end

function *(c::ZZRingElem, tdos::ToricDivisorOnScheme)
  result = ToricDivisorOnScheme(tdos.X, c*underlying_toric_divisor(tdos))
  if isdefined(tdos, :d)
    result.d = c*underlying_weil_divisor(tdos)
  end
  return result
end

function *(c::T, 
           tdos::ToricDivisorOnScheme{<:ToricCoveredScheme, <:AbstractAlgebra.Ring, T}
          ) where {T<:AbstractAlgebra.RingElem}
  result = ToricDivisorOnScheme(tdos.X, c*underlying_toric_divisor(tdos))
  if isdefined(tdos, :d)
    result.d = c*underlying_weil_divisor(tdos)
  end
  return result
end


function +(tdos1::ToricDivisorOnScheme, tdos2::ToricDivisorOnScheme)
  base_scheme(tdos1) === base_scheme(tdos2) || error("base schemes incompatible")
  return ToricDivisorOnScheme(base_scheme(tdos1), 
                              underlying_toric_divisor(tdos1) + underlying_toric_divisor(tdos2)
                             )
end

# In order to be able to also use a `ToricDivisorOnScheme` seamlessly 
# together with `ToricDivisor`s, we need to allow for the following 
# cross-type signatures.
function +(tdos1::ToricDivisorOnScheme, td2::ToricDivisor)
  toric_variety(tdos1) === toric_variety(td2) || error("toric varieties incompatible")
  return underlying_toric_divisor(tdos1) + td2
end
+(td2::ToricDivisor, tdos1::ToricDivisorOnScheme) = tdos1 + td2


function -(tdos1::ToricDivisorOnScheme, tdos2::ToricDivisorOnScheme)
  return tdos1 + (-1)*tdos2
end

function -(tdos1::ToricDivisorOnScheme, td2::ToricDivisor)
  return tdos1 + (-1)*td2
end

function -(td2::ToricDivisor, tdos1::ToricDivisorOnScheme)
  return td2 + (-1)*tdos1
end

function ==(tdos1::ToricDivisorOnScheme, tdos2::ToricDivisorOnScheme)
  base_scheme(tdos1) === base_scheme(tdos2) || error("base schemes incompatible")
  return underlying_toric_divisor(tdos1) == underlying_toric_divisor(tdos2)
end

function ==(tdos1::ToricDivisorOnScheme, td2::ToricDivisor)
  toric_variety(tdos1) === toric_variety(td2) || error("toric varieties incompatible")
  return underlying_toric_divisor(tdos1) == td2
end

function Base.show(io::IO, tdos::ToricDivisorOnScheme)
  return show(io, underlying_toric_divisor(tdos))
end

istrivial(tdos::ToricDivisorOnScheme) = istrivial(underlying_toric_divisor(tdos))
is_prime(tdos::ToricDivisorOnScheme) = is_prime(underlying_toric_divisor(tdos))
is_principal(tdos::ToricDivisorOnScheme) = is_principal(underlying_toric_divisor(tdos))
is_effective(tdos::ToricDivisorOnScheme) = is_effective(underlying_toric_divisor(tdos))
is_integral(tdos::ToricDivisorOnScheme) = is_integral(underlying_toric_divisor(tdos))
is_ample(tdos::ToricDivisorOnScheme) = is_ample(underlying_toric_divisor(tdos))
is_basepoint_free(tdos::ToricDivisorOnScheme) = is_basepoint_free(underlying_toric_divisor(tdos))
is_cartier(tdos::ToricDivisorOnScheme) = is_cartier(underlying_toric_divisor(tdos))
is_nef(tdos::ToricDivisorOnScheme) = is_nef(underlying_toric_divisor(tdos))
is_q_cartier(tdos::ToricDivisorOnScheme) = is_q_cartier(underlying_toric_divisor(tdos))
is_very_ample(tdos::ToricDivisorOnScheme) = is_very_ample(underlying_toric_divisor(tdos))

CohomologyClass(tdos::ToricDivisorOnScheme) = CohomologyClass(underlying_toric_divisor(tdos))
RationalEquivalenceClass(tdos::ToricDivisorOnScheme) = RationalEquivalenceClass(underlying_toric_divisor(tdos))
ToricDivisorClass(tdos::ToricDivisorOnScheme) = ToricDivisorClass(underlying_toric_divisor(tdos))
ToricLineBundle(tdos::ToricDivisorOnScheme) = ToricLineBundle(underlying_toric_divisor(tdos))

function ToricLineBundle(v::AbstractNormalToricVariety, tdos::ToricDivisorOnScheme) 
  return ToricLineBundle(toric_variety(underlying_toric_divisor(tdos)),
                         underlying_toric_divisor(tdos)
                        )
end

coefficients(tdos::ToricDivisorOnScheme) = coefficients(underlying_toric_divisor(tdos))
cohomology_class(tdos::ToricDivisorOnScheme) = cohomology_class(underlying_toric_divisor(tdos))
polyhedron(tdos::ToricDivisorOnScheme) = polyhedron(underlying_toric_divisor(tdos))
rational_equivalence_class(tdos::ToricDivisorOnScheme) = rational_equivalence_class(underlying_toric_divisor(tdos))
toric_divisor_class(tdos::ToricDivisorOnScheme) = toric_divisor_class(underlying_toric_divisor(tdos))
toric_line_bundle(tdos::ToricDivisorOnScheme) = toric_line_bundle(underlying_toric_divisor(tdos))
function toric_line_bundle(v::AbstractNormalToricVariety, tdos::ToricDivisorOnScheme) 
  return toric_line_bundle(v, underlying_toric_divisor(tdos))
end

toric_variety(tdos::ToricDivisorOnScheme) = toric_variety(underlying_toric_divisor(tdos))

