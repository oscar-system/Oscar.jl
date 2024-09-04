#function is_coxeter_matrix(M::ZZMatrix) end

@doc raw"""
    coxeter_from_cartan_matrix(mat::ZZMatrix; check::Bool=true) -> Bool

Return the Coxeter matrix $m$ associated to the Cartan matrix `gcm`. If there is no relation between $i$ and $j$,
then this will be expressed by $m_{ij} = 0$ (instead of the usual convention $m_{ij} = \infty$).
The keyword argument `check` can be set to `false` to skip verification whether `gcm` is indeed a generalized Cartan matrix.
"""
function coxeter_from_cartan_matrix(gcm; check::Bool=true)
  @req !check || is_cartan_matrix(gcm) "requires a generalized Cartan matrix"

  cm = identity_matrix(gcm)
  rk, _ = size(gcm)
  for i in 1:rk, j in 1:rk
    if i == j
      continue
    end

    if gcm[i, j] == 0
      cm[i, j] = 2
    elseif gcm[i, j] == -1
      cm[i, j] = 3
    elseif gcm[i, j] == -2
      cm[i, j] = 4
    elseif gcm[i, j] == -3
      cm[i, j] = 6
    end
  end

  return cm
end
