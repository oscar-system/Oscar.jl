###############################################################################
# A place to accumulate code that should eventually be moved to Nemo.jl
###############################################################################

function Hecke.numerator(f::QQPolyRingElem, parent::ZZPolyRing = Hecke.Globals.Zx)
  g = parent()
  ccall((:fmpq_poly_get_numerator, Nemo.libflint), Cvoid, (Ref{ZZPolyRingElem}, Ref{QQPolyRingElem}), g, f)
  return g
end

Hecke.minpoly(a::QQBarFieldElem) = minpoly(Hecke.Globals.Qx, a)

#TODO: Move to Nemo.
function minpoly(a::fpFieldElem)
  kx, x = polynomial_ring(parent(a), cached = false)
  return x-a
end
