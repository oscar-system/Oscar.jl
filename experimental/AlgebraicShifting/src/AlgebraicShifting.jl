function delta_ext(K::SimplicialComplex)
  n = nvertices(K)

  cmp(S1, S2) = min(symdiff(S1, S2)...) in S1
  K_facets = sort(facets(K); lt=cmp)
  m = length(K_facets)
  k = length(K_facets[1])
  
  K_min_non_faces = minimal_nonfaces(K)

  Qx, x = polynomial_ring(QQ, :x => (1:n, 1:n))
  F = fraction_field(Qx)
  X = matrix(F, hcat(x))
  nCk = sort(subsets(n, k))
  sub_compound_matrix = sparse_matrix(F, 0, 0)

  for row_subset in K_facets
    row_minors = Tuple{Int, AbstractAlgebra.Generic.FracFieldElem}[]

    for (i, col_subset) in enumerate(nCk)
      has_non_face = false
      for non_face in K_min_non_faces
        if is_subset(non_face, Set(col_subset))
          has_non_face = true
          break
        end
      end

      if has_non_face
        continue
      end
    
      ri = collect(row_subset)
      ci = collect(col_subset)
      push!(row_minors, (i, det(X[ri, ci])))
    end
    push!(sub_compound_matrix, sparse_row(F, row_minors))
  end

  @time _, _, _, U = lu(matrix(sub_compound_matrix));
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

