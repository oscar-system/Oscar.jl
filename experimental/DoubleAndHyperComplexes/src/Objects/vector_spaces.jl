# Some missing functionality to make complexes and double complexes of 
# vector spaces work.

function chain_complex(maps::Vector{T}; seed::Int=0, check::Bool=true) where {S<:FieldElem, T<:AbstractAlgebra.Generic.ModuleHomomorphism{S}}
  result = ComplexOfMorphisms(AbstractAlgebra.FPModule{S}, maps, typ=:chain, seed=seed, check=check)
  result.complete = true
  return result
end

function cochain_complex(maps::Vector{T}; seed::Int=0, check::Bool=true) where {S<:FieldElem, T<:AbstractAlgebra.Generic.ModuleHomomorphism{S}}
  result = ComplexOfMorphisms(AbstractAlgebra.FPModule{S}, maps, typ=:cochain, seed=seed, check=check)
  result.complete = true
  return result
end

morphism_type(::Type{T}) where {S<:FieldElem, T <: AbstractAlgebra.FPModule{S}} = AbstractAlgebra.Generic.ModuleHomomorphism{S}

zero_morphism(dom::AbstractAlgebra.FPModule, cod::AbstractAlgebra.FPModule) = hom(dom, cod, [zero(cod) for i in 1:ngens(dom)])

function Base.:*(k::Int, phi::AbstractAlgebra.Generic.ModuleHomomorphism)
  R = base_ring(codomain(phi))
  return R(k)*phi
end

# Tensor products of chain complexes of vector spaces
function (fac::TensorProductFactory{AbstractAlgebra.FPModule{S}})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) where {S<:FieldElem}
  return tensor_product(fac.C1[i], fac.C2[j])
end

# On input (dc, i, j) this produces the morphism 
#
#   φᵢ ⊗ id : Cᵢ⊗ Dⱼ → Cᵢ₊₋₁⊗ Dⱼ 
#
# induced by the (co-)boundary map φᵢ : Cᵢ → Cᵢ₊₋1 of the first complex 
# with the sign depending on the `horizontal_direction` of `dc`.
#
# See also experimental/src/tensor_products.jl
function (fac::HorizontalTensorMapFactory{AbstractAlgebra.Generic.ModuleHomomorphism{S}})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) where {S<:FieldElem}
  # construct the induced morphism C1[i] ⊗ C2[j] → C1[i ± 1] ⊗ C2[j]
  r2j = rank(fac.C2[j])
  A = matrix(map(fac.C1, i))

  # block_diagonal_matrix doesn't handle zero blocks
  if iszero(r2j) || iszero(ncols(A)) || iszero(nrows(A))
    inc = (typ(fac.C1) == :chain ? -1 : 1)
    return zero_morphism(dc[i, j], dc[i + inc, j])
  end

  res_mat = block_diagonal_matrix([A for i in 1:r2j])
  inc = (typ(fac.C1) == :chain ? -1 : 1)
  return hom(dc[i, j], dc[i + inc, j], res_mat)
end

# On input (dc, i, j) this produces the morphism 
#
#   id ⊗ ψⱼ : Cᵢ⊗ Dⱼ → Cᵢ⊗ Dⱼ₊₋₁
#
# induced by the (co-)boundary map ψⱼ : Dⱼ → Dⱼ₊₋1 of the second complex 
# with the sign depending on the `vertical_direction` of `dc`.
#
# See also experimental/src/tensor_products.jl
function (fac::VerticalTensorMapFactory{AbstractAlgebra.Generic.ModuleHomomorphism{S}})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int) where {S<:FieldElem}
  # construct the induced morphism C1[i] ⊗ C2[j] → C1[i] ⊗ C2[j ± 1]
  r1i = rank(fac.C1[i])
  B = matrix(map(fac.C2, j))
  
  # block_diagonal_matrix doesn't handle zero blocks
  if iszero(r1i) || iszero(ncols(B)) || iszero(nrows(B))
    inc = (typ(fac.C2) == :chain ? -1 : 1)
    return zero_morphism(dc[i, j], dc[i, j + inc])
  end

  R = base_ring(B)
  res_mat = zero_matrix(R, r1i*nrows(B), r1i*ncols(B))
  for ii in 1:nrows(B)
    for jj in 1:ncols(B)
      for k in 1:r1i
        res_mat[(ii - 1)*r1i + k, (jj - 1)*r1i + k] = B[ii, jj]
      end
    end
  end
  inc = (typ(fac.C2) == :chain ? -1 : 1)
  return hom(dc[i, j], dc[i, j + inc], res_mat)
end


