function w_to_alpha(
  L::LieAlgebraStructure, weight_w::Union{Vector{ZZRingElem},Vector{QQFieldElem}}
)
  return weight_w * inv_cartan_matrix(L)
end

function alpha_to_w(L::LieAlgebraStructure, weight_alpha::Vector{QQFieldElem})
  return weight_alpha * cartan_matrix(L)
end
