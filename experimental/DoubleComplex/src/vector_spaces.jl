# Some missing functionality to make complexes and double complexes of 
# vector spaces work.

function chain_complex(maps::Vector{T}; seed::Int=0, check::Bool=true) where {S<:FieldElem, T<:AbstractAlgebra.Generic.ModuleHomomorphism{S}}
  result = ComplexOfMorphisms(AbstractAlgebra.Generic.FreeModule{S}, maps, typ=:chain, seed=seed, check=check)
  result.complete = true
  return result
end

function cochain_complex(maps::Vector{T}; seed::Int=0, check::Bool=true) where {S<:FieldElem, T<:AbstractAlgebra.Generic.ModuleHomomorphism{S}}
  result = ComplexOfMorphisms(AbstractAlgebra.Generic.FreeModule{S}, maps, typ=:cochain, seed=seed, check=check)
  result.complete = true
  return result
end

Oscar.morphism_type(::Type{T}) where {S<:FieldElem, T <: AbstractAlgebra.Generic.FreeModule{S}} = AbstractAlgebra.Generic.ModuleHomomorphism{S}

zero_morphism(dom::AbstractAlgebra.Generic.FreeModule, cod::AbstractAlgebra.Generic.FreeModule) = hom(dom, cod, [zero(cod) for i in 1:ngens(dom)])

# Tensor products of chain complexes of vector spaces
function (fac::Oscar.TensorProductFactory{AbstractAlgebra.Generic.FreeModule{S}})(dc::Oscar.AbsDoubleComplexOfMorphisms, i::Int, j::Int) where {S<:FieldElem}
  return tensor_product(fac.C1[i], fac.C2[j])
end

function (fac::Oscar.HorizontalTensorMapFactory{AbstractAlgebra.Generic.ModuleHomomorphism{S}})(dc::Oscar.AbsDoubleComplexOfMorphisms, i::Int, j::Int) where {S<:FieldElem}
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

function (fac::Oscar.VerticalTensorMapFactory{AbstractAlgebra.Generic.ModuleHomomorphism{S}})(dc::Oscar.AbsDoubleComplexOfMorphisms, i::Int, j::Int) where {S<:FieldElem}
  # construct the induced morphism C1[i] ⊗ C2[j] → C1[i] ⊗ C2[j ± 1]
  r1i = rank(fac.C1[i])
  B = matrix(map(fac.C2, j))
  
  # block_diagonal_matrix doesn't handle zero blocks
  if iszero(r1i) || iszero(ncols(B)) || iszero(nrows(B))
    inc = (typ(fac.C2) == :chain ? -1 : 1)
    return zero_morphism(dc[i, j], dc[i, j + inc])
  end

  R = base_ring(B)
  res_mat = zero(MatrixSpace(R, r1i*nrows(B), r1i*ncols(B)))
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

