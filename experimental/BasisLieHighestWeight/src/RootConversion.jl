function w_to_alpha(
  L::LieAlgebraStructure, weight_w::Union{Vector{ZZRingElem},Vector{QQFieldElem}}
)
  return cartan_matrix_inv(L) * weight_w
end

function alpha_to_w(L::LieAlgebraStructure, weight_alpha::Vector{QQFieldElem})
  return cartan_matrix(L) * weight_alpha
end
