#function is_coxeter_matrix(M::ZZMatrix) end

function coxeter_matrix(W::WeylGroup)
  return cartan_to_coxeter_matrix(cartan_matrix(W))
end

@doc raw"""
    coxeter_from_cartan_matrix(mat::ZZMatrix; check::Bool=true) -> Bool

Return the Coxeter matrix $m$ associated to the Cartan matrix `gcm`. If there is no relation between $i$ and $j$,
then this will be expressed by $m_{ij} = 0$ (instead of the usual convention $m_{ij} = \infty$).
The keyword argument `check` can be set to `false` to skip verification whether `gcm` is indeed a generalized Cartan matrix.
"""
function coxeter_from_cartan_matrix(gcm::ZZMatrix; check::Bool=true)
  @req !check || is_cartan_matrix(gcm) "requires a generalized Cartan matrix"

  rk, _ = size(gcm)
  cm = matrix(
    ZZ, [coxeter_matrix_entry_from_cartan_matrix(gcm, i, j) for i in 1:rk, j in 1:rk]
  )

  return cm
end

function coxeter_matrix_entry_from_cartan_matrix(gcm::ZZMatrix, i::Int, j::Int)
  if i == j
    return 1
  end

  d = gcm[i, j] * gcm[j, i]
  if d == 0
    return 2
  elseif d == 1
    return 3
  elseif d == 2
    return 4
  elseif d == 3
    return 6
  else
    return 0
  end
end
