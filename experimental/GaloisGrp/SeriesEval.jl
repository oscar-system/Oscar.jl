function (f::RelPowerSeriesRingElem)(a::RingElem)
  y = a
  v = valuation(f)
  p = precision(f)
  z = zero(y)
  for i = Int(p):-1:Int(v)
      z *= y
      c = coeff(f, i)
      if !iszero(c)
          z += c
      end
  end
  z *= y^Int(v)
  return z
end

