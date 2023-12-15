function delta_ext(K::SimplicialComplex)
  n = nvertices(K)

  cmp(S1, S2) = min(symdiff(S1, S2)...) in S1
  K_facets = sort(facets(K); lt=cmp)

  k = length(K_facets[1])
  Qx, x = polynomial_ring(QQ, :x => (1:n, 1:n))
  X = matrix(fraction_field(Qx), hcat(x))
  sub_compound_matrix_rows = Vector{FracElem{<:MPolyRingElem}}[]
  nCk = sort(subsets(n, k))
  
  for col_subset in nCk
    col_minors = []
    for row_subset in K_facets
      ri = collect(row_subset)
      ci = collect(col_subset)
      push!(col_minors, det(X[ri, ci]))
    end
    push!(sub_compound_matrix_rows, col_minors)
  end

  A = reduce(hcat, sub_compound_matrix_rows)
  _, _, _, U = lu(matrix(A));
  delta_K = Int[]

  row_index = 1
  row_bound = nrows(A)
  for col_index in 1:ncols(A)
    if !is_zero(U[row_index, col_index])
      push!(delta_K, col_index)
      row_index += 1
    end

    row_index > row_bound && break
  end
  return simplicial_complex(nCk[delta_K])
end

#@doc raw"""
#    my_access_func(S::ExampleStruct)
#
#This is a dummy sample function to teach the use of docstrings.
#"""
#function my_access_func(S::ExampleStruct)
#  return S.i
#end

