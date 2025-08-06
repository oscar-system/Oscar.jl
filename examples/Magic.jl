module Magic

using Oscar

function magic_square_generators(n::Int)
  l = Vector{Vector{Int}}()
  for i=1:1
    for j=i+1:n
      m = zeros(Int, n*n)
      #sum row[1] == sum row[i]
      for k=1:n
        m[(i-1)*n+k] = 1
        m[(j-1)*n+k] = -1
      end
      push!(l, m)
    end
    #sum row[1] = sum col[i] (so [1,i] is zero: plus 1, minus 1)
    for j=1:n
      m = zeros(Int, n*n)
      for k=1:n
        m[(i-1)*n+k] += 1
        m[(k-1)*n+j] -= 1
      end
      push!(l, m)
    end
    #sum row[1]= sum diag
    m = zeros(Int, n*n)
    for k=2:n
      m[k] = 1
      m[(k-1)*n+k] = -1
    end
    push!(l, m)
    #sum row[1]= sum anti-diag
    m = zeros(Int, n*n)
    for k=1:n-1
      m[k] = 1
      m[(n-k)*n+k] = -1
    end
    push!(l, m)
  end
  l = matrix(ZZ, l)
  P = Oscar._cone(A = l, b = zero_matrix(ZZ, nrows(l), 1), C = identity_matrix(ZZ, ncols(l)))
  HB = P.HILBERT_BASIS_GENERATORS
  return matrix(ZZ, HB[1])
  #follows https://www.researchgate.net/publication/2101091_Polyhedral_Cones_of_Magic_Cubes_and_Squares
  #next step requires quotient rings (and possibly subquotient rings)
end

function Oscar.matrix(::ZZRing, M::Polymake.MatrixAllocated{Polymake.Integer})
  s = size(M)
  N = zero_matrix(ZZ, s[1], s[2]-1)
  for i=1:s[1]
    for j=2:s[2]
      N[i, j-1] = ZZRingElem(M[i,j])
    end
  end
  return N
end

end
