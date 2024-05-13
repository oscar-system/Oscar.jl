
function _to_integer_column_matrix(v::SubObjectIterator{PointVector{ZZRingElem}})
  n = length(v)
  @assert n > 0 "can not convert the empty list"
  m = length(first(v))
  result = zero_matrix(ZZ, m, n)
  for j in 1:n
    for i in 1:m
      result[i, j] = v[j][i]
    end
  end
  return result
end

function _colum_vectors_to_rays(K::QQMatrix)
  m = nrows(K)
  n = ncols(K)
  result = Vector{Vector{QQFieldElem}}()
  for i in 1:n
    row = Vector{QQFieldElem}()
    for j in 1:m
      push!(row, K[j, i])
    end
    push!(result, row)
  end
  return result
end

function _colum_vectors_to_rays(K::ZZMatrix)
  m = nrows(K)
  n = ncols(K)
  result = Vector{Vector{ZZRingElem}}()
  for i in 1:n
    row = Vector{ZZRingElem}()
    for j in 1:m
      push!(row, K[j, i])
    end
    push!(result, row)
  end
  return result
end



function _to_matrix(v::RayVector{QQFieldElem})
  n = length(v)
  result = zero_matrix(QQ, 1, n)
  for i in 1:n
    result[1, i] = v[i]
  end
  return result
end

function _to_matrix(v::SubObjectIterator{RayVector{QQFieldElem}})
  m = length(v)
  @assert m > 0 "can not convert the empty list"
  n = length(first(v))
  @assert all(x->length(x)==n, v) "vectors must have the same lengths"
  result = zero_matrix(QQ, m, n)
  for i in 1:m
    for j in 1:n
      result[i, j] = v[i][j]
    end
  end
  return result
end