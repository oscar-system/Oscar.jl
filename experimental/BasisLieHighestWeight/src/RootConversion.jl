function w_to_alpha(type::Symbol, rank::Int, weight_w::Vector{QQFieldElem})
  return weight_w * inv(cartan_matrix(type, rank))
end

function alpha_to_w(type::Symbol, rank::Int, weight_alpha::Vector{QQFieldElem})
  return weight_alpha * cartan_matrix(type, rank)
end

function cartan_matrix(type::Symbol, rank::Int)
  L = GAP.Globals.SimpleLieAlgebra(GAP.Obj(type), rank, GAP.Globals.Rationals)
  R = GAP.Globals.RootSystem(L)
  C = matrix(QQ, GAP.Globals.CartanMatrix(R))
  return C
end
