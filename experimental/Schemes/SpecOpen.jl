export adjoint

#(R::MPolyQuoLocalizedRing)(f::SpecOpenRingElem) = restrict(f, Spec(R))
#(R::MPolyLocalizedRing)(f::SpecOpenRingElem) = restrict(f, Spec(R))
#(R::MPolyRing)(f::SpecOpenRingElem) = restrict(f, Spec(R))


#function Base.adjoint(M::AbstractAlgebra.Generic.MatSpaceElem{T}) where {T} 
function Base.adjoint(M::MatElem)
  n = ncols(M)
  n == nrows(M) || error("matrix is not square")
  N = zero(M)
  rows = collect(1:n)
  cols = collect(1:n)
  row_sign = 1
  for i in 1:n
    column_sign = row_sign
    for j in 1:n
      N[j, i] = column_sign* det(M[deleteat!(copy(rows), i), deleteat!(copy(cols), j)])
      column_sign = column_sign*(-1)
    end
    row_sign = row_sign*(-1)
  end
  return N
end

Base.inv(M::MatElem) = inv(det(M))*adjoint(M)


  
