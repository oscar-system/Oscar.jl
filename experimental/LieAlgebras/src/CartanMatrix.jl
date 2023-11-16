###############################################################################
#
#   Cartan Matrix Helpers
#   
###############################################################################

@doc raw"""
    cartan_matrix(fam::Symbol, rk::Int) -> ZZMatrix
"""
function cartan_matrix(fam::Symbol, rk::Int)
  if fam == :A
    @req rk >= 1 "type An requires rank rk of at least 1"

    mat = diagonal_matrix(ZZ(2), rk)
    for i in 1:(rk - 1)
      mat[i + 1, i], mat[i, i + 1] = -1, -1
    end
  elseif fam == :B
    @req rk >= 2 "type Bn requires rank rk of at least  2"

    mat = diagonal_matrix(ZZ(2), rk)
    for i in 1:(rk - 1)
      mat[i + 1, i], mat[i, i + 1] = -1, -1
    end
    mat[rk, rk - 1] = -2
  elseif fam == :C
    @req rk >= 2 "type Cn requires rank rk of at least 2"

    mat = diagonal_matrix(ZZ(2), rk)
    for i in 1:(rk - 1)
      mat[i + 1, i], mat[i, i + 1] = -1, -1
    end
    mat[rk - 1, rk] = -2
  elseif fam == :D
    @req rk >= 3 "type Dn requires rank rk of at least 3"

    mat = diagonal_matrix(ZZ(2), rk)
    for i in 1:(rk - 1)
      mat[i + 1, i], mat[i, i + 1] = -1, -1
    end
    mat[rk - 1, rk] = -2
  elseif fam == :E
    @req rk in 6:8 "type En requires rank to be one of 6, 7, 8"
    mat = matrix(
      ZZ,
      [
        2 0 -1 0 0 0 0 0
        0 2 0 -1 0 0 0 0
        -1 0 2 -1 0 0 0 0
        0 -1 -1 2 -1 0 0 0
        0 0 0 -1 2 -1 0 0
        0 0 0 0 -1 2 -1 0
        0 0 0 0 0 -1 2 -1
        0 0 0 0 0 0 -1 2
      ],
    )

    if rk == 6
      mat = mat[1:6, 1:6]
    elseif rk == 7
      mat = mat[1:7, 1:7]
    end
  elseif fam == :F
    @req rk == 4 "type Fn requires rank rk to be 4"
    mat = matrix(ZZ, [2 -1 0 0; -1 2 -1 0; 0 -2 2 -1; 0 0 -1 2])
  elseif fam == :G
    @req rk == 2 "type Gn requires rank rk to be 2"
    mat = matrix(ZZ, [2 -3; -1 2])
  else
    error("unknown family $fam")
  end

  return mat
end

function cartan_matrix(types::Tuple{Symbol,Int}...)
  blocks = ZZMatrix[cartan_matrix(type...) for type in types]

  return block_diagonal_matrix(blocks)
end

function cartan_to_coxeter_matrix(mat; check=true)
  @req !check || !is_cartan_matrix(mat) "$mat is not a (generalized) Cartan matrix"

  cm = identity_matrix(mat)
  rk, _ = size(mat)
  for i in 1:rk, j in 1:rk
    if i == j
      continue
    end

    if mat[i, j] == 0
      cm[i, j] = 2
    elseif mat[i, j] == -1
      cm[i, j] = 3
    elseif mat[i, j] == -2
      cm[i, j] = 4
    elseif mat[i, j] == -3
      cm[i, j] = 6
    end
  end

  return cm
end

@doc raw"""
    is_cartan_matrix(mat::ZZMatrix) -> Bool

Test if `mat` is a (`generalized`) Cartan matrix.
"""
function is_cartan_matrix(mat::ZZMatrix; generalized::Bool=true)
  n, m = size(mat)
  if n != m
    return false
  end

  if !all(mat[i, i] == 2 for i in 1:n)
    return false
  end

  for i in 1:n
    for j in (i + 1):n
      if is_zero_entry(mat, i, j) && is_zero_entry(mat, j, i)
        continue
      end

      # mat[i,j] != 0 or mat[j,i] != 0, so both entries must be < 0
      if mat[i, j] >= 0 || mat[j, i] >= 0
        return false
      end

      # if we only want a generalized Cartan matrix, we are done with this pair of entries
      if generalized
        continue
      end

      if !(mat[i, j] * mat[j, i] in 1:3)
        return false
      end
    end
  end

  return true
end
