
function produce_object_on_affine_chart(II::ToricIdealSheafFromCoxRingIdeal, U::AbsAffineScheme)
  X = scheme(II)
  k = findfirst(V->V===U, affine_charts(X))
  k === nothing && error("affine chart not found")
  return _dehomogenize_to_chart(X, ideal_in_cox_ring(II), k)
end
