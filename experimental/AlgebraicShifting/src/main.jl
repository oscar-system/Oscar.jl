# Add your new types, functions, and methods here.
n = 5
k = 3
m = 4

Qx, x = polynomial_ring(QQ, :x => (1:n, 1:n))

X = matrix(fraction_field(Qx), hcat(x))

possible_m_sets = subsets(binomial(n, k), m)

for i in 1:length(possible_m_sets)
  K = sort(subsets(n, k)[possible_m_sets[i]])

  sub_compound_matrix_rows = Vector{AbstractAlgebra.Generic.Frac}[]
  for col_subset in sort(subsets(n, k))
    col_minors = []
    for row_subset in K
      push!(col_minors, det(X[row_subset, col_subset]))
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

  println(delta_K)
end

@doc raw"""
    my_access_func(S::ExampleStruct)

This is a dummy sample function to teach the use of docstrings.
"""
function my_access_func(S::ExampleStruct)
  return S.i
end

