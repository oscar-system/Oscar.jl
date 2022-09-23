export AffineMatrixGroup
export coordinate_matrix, coordinate_ring, modulus, inverted_set, ground_field
export special_linear_group

@attributes mutable struct AffineMatrixGroup{BaseRingType, RingType, 
                                       SpecType<:AbsAffineGroupScheme, 
                                       MatrixType<:MatrixElem,
                                       FuncType1, FuncType2
                                      } <: AbsAffineGroupScheme{BaseRingType, RingType}
  G::SpecType # the underlying group scheme
  M::MatrixType# the coordinate matrix
  to_linear_index::FuncType1 # translation between the indices of the coordinates and the matrix entries
  to_matrix_index::FuncType2
  id_ideal::MPolyIdeal
  
  function AffineMatrixGroup(M::MatrixElem{T}, I::MPolyIdeal{T}; 
      check::Bool=true
    ) where {T<:MPolyElem}

    R = base_ring(M)
    m = ncols(M)
    m == nrows(M) || error("coordinate matrix is not square")
    R == base_ring(I) || error("matrix and ideal are not defined over the same ring")
    kk = coefficient_ring(R)

    # check that the coordinate matrix has only coordinates as entries
    # and assemble the transition functions for the indices
    a = copy(gens(R))
    index_lin = Array{Int}(undef, m, m)
    index_mat = [[0, 0] for i in 1:ngens(R)]
    for i in 1:m
      for j in 1:m
        tmp = findfirst(x->(isequal(x, M[i,j])), a)
        tmp == nothing && error("the coordinate matrix must have pairwise different coordinate functions as entries")
        index_lin[i, j] = tmp
        index_mat[tmp] = Int[i, j]
        a[tmp] = zero(R)
      end
    end
    to_linear(i::Int, j::Int) = index_lin[i,j]
    to_square(k::Int) = index_mat[k]

    # Set up the various maps needed to create a group scheme.
    # The product of X with itself:
    f = det(M)
    X = Spec(R, I, MPolyPowersOfElement(R, [f]))
    XxX, p1, p2 = product(X, X, change_var_names_to=["y", "z"])

    # The diagonal embedding:
    pb_diag = hom(OO(XxX), OO(X), vcat(gens(R), gens(R)))
    diag = SpecMor(X, XxX, pb_diag)

    # The multiplication map:
    A = map_entries(pullback(p1), M)
    B = map_entries(pullback(p2), M)
    C = A*B
    pb_prod_map = hom(OO(X), OO(XxX), [C[to_square(k)[1], to_square(k)[2]] for k in 1:ngens(R)])
    prod_map = SpecMor(XxX, X, pb_prod_map, check=check)

    # The vanishing ideal of the unit element for caching:
    id_ideal = ideal(R, [gens(R)[k] - (to_square(k)[1] == to_square(k)[2] ? 1 : 0) for k in 1:ngens(R)])

    # The inclusion maps for the neutral element in the other factor:
    pb_inc1 = hom(OO(XxX), OO(X), 
                  vcat(gens(R), R.([(to_square(k)[1] == to_square(k)[2] ? 1 : 0) for k in 1:ngens(R)])))
    inc1 = SpecMor(X, XxX, pb_inc1, check=check)
    pb_inc2 = hom(OO(XxX), OO(X), 
                  vcat(R.([(to_square(k)[1] == to_square(k)[2] ? 1 : 0) for k in 1:ngens(R)]), gens(R)))
    inc2 = SpecMor(X, XxX, pb_inc2, check=check)

    # The inversion map for the group law:
    inv_M = _help_inv(map_entries(OO(X), M))
    pb_inv = hom(OO(X), OO(X), [inv_M[to_square(k)[1], to_square(k)[2]] for k in 1:ngens(R)])
    inv_map = SpecMor(X, X, pb_inv, check=check)

    G = AffineGroupScheme(X, XxX, diag, p1, p2, inc1, inc2, prod_map, inv_map, 
                          kk.([(to_square(k)[1] == to_square(k)[2] ? 1 : 0) for k in 1:ngens(R)]),
                          check=check)
    return new{typeof(kk), typeof(OO(G)), 
               typeof(G), typeof(M), 
               typeof(to_linear), typeof(to_square)
              }(
                G, M, to_linear, to_square, id_ideal
               )

   
  end
end

### essential getters
underlying_group_scheme(G::AffineMatrixGroup) = G.G
underlying_scheme(G::AffineMatrixGroup) = underlying_scheme(G.G)

### user facing essential getters
coordinate_matrix(G::AffineMatrixGroup) = G.M
to_linear_index_func(G::AffineMatrixGroup) = G.to_linear_index
to_matrix_index_func(G::AffineMatrixGroup) = G.to_matrix_index
ideal_of_neutral_element(G::AffineMatrixGroup) = G.id_ideal

### derived getters
coordinate_ring(G::AffineMatrixGroup) = base_ring(coordinate_matrix(G))
modulus(G::AffineMatrixGroup) = modulus(OO(G))
inverted_set(G::AffineMatrixGroup) = inverted_set(OO(G))
ground_field(G::AffineMatrixGroup) = coefficient_ring(base_ring(coordinate_matrix(G)))

### 
# Given a matrix A and an affine matrix group G, this function 
# returns the coordinate vector for the point A in the coordinates 
# of G (the variables of the underlying polynomial ring).
function coordinates(G::AffineMatrixGroup, A::MatrixElem)
  M = coordinate_matrix(G)
  m = ncols(M) 
  m == ncols(A) == nrows(A) || error("matrix sizes are not compatible")
  ind = to_matrix_index_func(G)
  return [evaluate(M[ind(k)[1], ind(k)[2]], A[ind(k)[1], ind(k)[2]]) 
          for k in 1:ngens(coordinate_ring(G))]
end
                   

function _help_inv(M::MatrixElem)
  m = ncols(M)
  m == nrows(M) || error("inverse can only be computed for square matrices")
  N = zero(M)
  row_sign = 1
  for i in 1:m
    sign = row_sign
    for j in 1:m
      N[i,j] = sign*det(
                        vcat(hcat(M[1:i-1, 1:j-1], M[1:i-1, j+1:m]),
                             hcat(M[i+1:m, 1:j-1], M[i+1:m, j+1:m]))
                       )
      sign = -sign
    end
    row_sign = - row_sign
  end
  f = det(M)
  f = inv(f)
  return f*N
end

function special_linear_group(::Type{AffineMatrixGroup}, n::Int, kk::Field)
  R, a = PolynomialRing(kk, ["a_$(i)_$(j)" for i in 1:n for j in 1:n])
  M = zero(MatrixSpace(R, n, n))
  for i in 1:n
    for j in 1:n
      M[i,j] = a[(i-1)*n+j]
    end
  end
  f = det(M)
  I = ideal(R, [f-1])
  return AffineMatrixGroup(M, I, check=true) #set to false, eventually!!!
end

function general_linear_group(::Type{AffineMatrixGroup}, n::Int, kk::Field)
  R, a = PolynomialRing(kk, ["a_$(i)_$(j)" for i in 1:n for j in 1:n])
  M = zero(MatrixSpace(R, n, n))
  for i in 1:n
    for j in 1:n
      M[i,j] = a[(i-1)*n+j]
    end
  end
  f = det(M)
  return AffineMatrixGroup(M, ideal(R, zero(R)), check=true) #set to false, eventually!!!
end


