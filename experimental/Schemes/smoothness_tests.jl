import AbstractAlgebra.Generic: MatSpaceElem

export degeneracy_locus, is_equidimensional_and_smooth, singular_locus

function degeneracy_locus(
    X::SpecType,
    A::MatrixType,
    r::Int; 
    rec_count::Int=0,
    verbose::Bool=false
  ) where {SpecType<:Spec, MatrixType} #TODO: How specific can we be with the matrices?
  indent_str = prod([ "#" for i in 1:rec_count ]) * " "
  verbose && println(indent_str * "call with $(X) and matrix $(A) for rank $r")

  # checking for rank < 0 does not make sense unless we have the the zero matrix
  if r == 0
    for a in collect(A)
      iszero(OO(X)(a)) || error("input not admissible")
    end
    return SpecType[]
  end

  # if X is already empty, quit
  isempty(X) && return SpecType[]

  # check for some other aborting conditions
  m = nrows(A)
  n = ncols(A)
  R = base_ring(OO(X))
  if r == 1 
    # if the matrix is trivial, but the rank is not, quit.
    if m == 0 || n == 0
      verbose && println(indent_str, "zero matrix; returning $X")
      return [X]
    end
    # if we are in rank < 1 for a non-trivial matrix, collect the entries and quit
    Y = subscheme(X, ideal(R, [A[i,j] for i in 1:nrows(A) for j in 1:ncols(A)]))
    isempty(Y) && return SpecType[]
    return [Y]
  end
  
  # find the entry with lowest degree to cut out the next hypersurface
  (k, l) = (1, 1)
  d = maximum(total_degree.(A))
  I = localized_modulus(OO(X))
  allzero = true
  W = localized_ring(OO(X))
  for i in 1:m
    verbose && (print_str = indent_str)
    for j in 1:n
      A[i,j] = numerator(reduce(W(A[i,j]), groebner_basis(I))) # TODO: Implement and use reduction for matrices?
      verbose && (print_str *= "$(A[i, j]);\t")
      allzero = allzero && iszero(A[i,j])
      if total_degree(A[i,j]) <= d && !iszero(A[i,j])
	d = total_degree(A[i,j])
	(k, l) = (i, j)
      end
    end
    verbose && println(print_str)
  end
  f = A[k,l]

  # in case that after reduction all the matrix' entries are zero, quit
  if allzero
    verbose && println(indent_str * "All entries are zero; returning $X")
    isempty(X) && return SpecType[]
    return [X]
  end
  verbose && println(indent_str * "selected entry at ($k, $l): $f")
  U = hypersurface_complement(X, f)
  Y = subscheme(X, f)
  B = copy(A)
  for i in 1:k-1
    multiply_row!(B, f, i)
    add_row!(B, -A[i,l], k, i)
  end
  for i in k+1:m
    multiply_row!(B, f, i)
    add_row!(B, -A[i,l], k, i)
  end

  # set up the submatrix for induction on the open part
  u = [i for i in 1:m if i != k]
  v = [j for j in 1:n if j != l]
  B = A[u, v]
  #B = vcat(hcat(B[1:k-1, 1:l-1], B[1:k-1, l+1:n]), hcat(B[k+1:m, 1:l-1], B[k+1:m, l+1:n]))
  #
  verbose && println(indent_str * "new matrix of size $(nrows(B)) x $(ncols(B)):")
  if verbose
    for i in 1:nrows(B)
      print_str = indent_str
      for j in 1:ncols(B)
        print_str = print_str * "$(B[i,j]);\t"
      end
      println(print_str)
    end
  end
  
  DU = degeneracy_locus(U, B, r-1, rec_count=rec_count+1)
  DU_closure = SpecType[closure(V, X) for V in DU]
  DY = degeneracy_locus(Y, copy(A), r, rec_count=rec_count+1)
  return vcat(DU_closure, DY)
end

function is_equidimensional_and_smooth(X::Spec; verbose::Bool=false)
  R = base_ring(OO(X))
  W = localized_ring(OO(X))
  IW = localized_modulus(OO(X))
  I = saturated_ideal(IW)
  d = dim(I)
  n = nvars(R)
  return length(degeneracy_locus(X, jacobi_matrix(gens(I)), n-d, verbose=verbose))== 0
end

function singular_locus(X::Spec; verbose::Bool=false)
  R = base_ring(OO(X))
  W = localized_ring(OO(X))
  IW = localized_modulus(OO(X))
  I = saturated_ideal(IW)
  d = dim(I)
  n = nvars(R)
  list = degeneracy_locus(X, jacobi_matrix(gens(I)), n-d, verbose=verbose)
  if length(list) == 0
    return hypersurface_complement(X, one(base_ring(OO(X))))
  else
    return subscheme(X, intersect(saturated_ideal.(localized_modulus.(OO.(list)))...))
  end
end
  

#TODO: This is not efficient yet, because of double checks!
function is_equidimensional_and_smooth(X::CoveredScheme)
  C = default_covering(X) 
  for U in patches(C)
    is_equidimensional_and_smooth(U) || return false
  end
  return true
end
