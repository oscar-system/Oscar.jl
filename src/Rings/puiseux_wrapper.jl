function _puiseux(f::MPolyRingElem{T}, max_deg::Int, at_origin::Bool=true) where {T<:FieldElem}
  P = parent(f)
  SP = singular_poly_ring(P)
  Sf = SP(f)
  return Singular.LibPuiseuxexpansions.puiseux(Sf, max_deg, Int(at_origin))
end
