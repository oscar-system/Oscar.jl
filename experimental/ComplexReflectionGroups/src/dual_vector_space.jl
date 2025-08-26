# Some functions concerning the dual of a vector space.
#
# Ulrich Thiel, 2023 


function canonical_pairing(v::AbstractAlgebra.Generic.FreeModuleElem{T}, w::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: FieldElem

  V = parent(v)
  K = base_ring(V)
  n = dim(V)
  s = zero(K)
  for i=1:n
    s += v[i]*w[i]
  end
  return s
end

function linear_form(v::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: FieldElem
   
  V = parent(v)
  K = base_ring(V)
  n = dim(V)

  form_matrix = matrix(K,n,1,[v[i] for i=1:n])
  form_hom = module_homomorphism(V, vector_space(K,1), form_matrix)

  return form_hom

end
