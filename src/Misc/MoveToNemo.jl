###############################################################################
# A place to accumulate code that should eventually be moved to Nemo.jl
###############################################################################

function minpoly(a::fpFieldElem)
  kx, x = polynomial_ring(parent(a), cached = false)
  return x-a
end
