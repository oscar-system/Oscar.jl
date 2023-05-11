mutable struct ToricDivisorOnScheme{CoveredSchemeType<:ToricCoveredScheme,
                                    CoeffRingType<:AbstractAlgebra.Ring, 
                                    CoeffRingElemType<:AbstractAlgebra.RingElem
                                   }<:AbsAlgebraicCycle{CoveredSchemeType, CoeffRingType}

  td::ToricDivisor
  X::CoveredSchemeType
  d::WeilDivisor{CoveredSchemeType, CoeffRingType, CoeffRingElemType}

  function ToricDivisorOnScheme(X::ToricCoveredScheme, td::ToricDivisor)
    return new{typeof(X), typeof(ZZ), elem_type(ZZ)}(td, X)
  end

  function ToricDivisorOnScheme(td::ToricDivisor)
    return new{ToricCoveredScheme, typeof(ZZ), elem_type(ZZ)}(td)
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

function base_scheme(tdos::ToricDivisorOnScheme)
  if !isdefined(tdos, :X)
    X = ToricCoveredScheme(toric_variety(underlying_toric_divisor(tdos)))
    tdos.X = X
  end
  return tdos.X
end

### inherited functionality from the toric divisors 
function *(c::T, tdos::ToricDivisorOnScheme) where {T<:IntegerUnion}
  if isdefined(tdos, :X)
    # If some data has already been computed, keep it.
    result = ToricDivisorOnScheme(tdos.X, c*underlying_toric_divisor(tdos))
    if isdefined(tdos, :d)
      result.d = c*underlying_weil_divisor(tdos)
    end
    return result
  else
    # If the base scheme is not set, then there can also be no Weil divisor, 
    # as this would have done the job.
    return ToricDivisorOnScheme(c*underlying_toric_divisor(tdos))
  end
end

# More to come!
