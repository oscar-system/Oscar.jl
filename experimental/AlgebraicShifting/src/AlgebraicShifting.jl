@doc raw"""
   delta_ext(K::SimplicialComplex)

Returns the exterior shift of `K`
"""
function delta_ext(K::SimplicialComplex)
  n = nvertices(K)
  dim_K = dim(K)

  cmp(S1, S2) = min(symdiff(S1, S2)...) in S1
  K_facets = sort(facets(K); lt=cmp)

  facet_dict = Dict(
    (k => filter(f -> length(f) == k + 1, K_facets) for k in 1:dim_K )
  )

  for (dim_face, faces) in facet_dict
    Qx, x = polynomial_ring(QQ, :x => (1:n, 1:n))
    X = matrix(Qx, hcat(x))
    nCk = sort(subsets(n, dim_face + 1))
    sub_compound_matrix = Vector{MPolyRingElem}[]

    for col_subset in nCk
      row_minors = MPolyRingElem[]
      for row_subset in faces
        ri = collect(row_subset)
        ci = collect(col_subset)
        push!(row_minors, det(X[ri, ci]))
      end
      push!(sub_compound_matrix, row_minors)
    end

    A = matrix(Qx, sub_compound_matrix)
    Oscar.ModStdQt.ref_ff_rc!(A)

    delta_k = Int[]
    row_index = 1
    row_bound = nrows(A)
    for col_index in 1:ncols(A)
      if !is_zero(A[row_index, col_index])
        push!(delta_k, col_index)
        row_index += 1
      end

      row_index > row_bound && break
    end
    facet_dict[dim_face] = map(Set, nCk[delta_k])
  end
  return simplicial_complex(
    reduce(vcat, values(facet_dict))
  )
end

@doc raw"""
   delta_sym(K::SimplicialComplex)

Returns the symmetric shift of `K`
"""
function delta_sym(K::SimplicialComplex)
  R_K, quotient_map = stanley_reisner_ring(K)
  Ry, y = polynomial_ring(R_K, :y => (1:n, 1:n))
  Y = matrix(Ry, hcat(y))
  y_1 = Y[1, :] * gens(R_K)
  
  collect(terms(y_1[1] * y_1[1]))
end

function generic_initial_ideal(I::MPolyIdeal)
  R = base_ring(I)
  n = nvars(R)
  Ry, y = polynomial_ring(R, :y => (1:n, 1:n))
  Y = matrix(Ry, hcat(y))
  gens_R = gens(R)
  generic_ini_gens = MPolyRingElem[]

  for g in gens(I)
    monomials = exponents(g)
    non_zero_exp_indices = (findall(x -> !is_zero(x), m) for m in monomials)
    generic_gen = zero(Ry)
    for (m, e) in zip(monomials, non_zero_exp_indices)
      degrees = m[e]
      multinomials = Y[e, :] * gens_R
      generic_monomial = prod(zip(multinomials, degrees)) do (multinomial, d)
        multinomial^d
      end
      generic_gen += generic_monomial
    end
    push!(generic_ini_gens, generic_gen)
  end
  return ideal(generic_ini_gens)
end
