include("UniformHypergraph.jl")
include("PartialShift.jl")
include("PartialShiftGraph.jl")

export UniformHypergraph

export compound_matrix
export contracted_partial_shift_graph
export exterior_shift
export face_size
export generic_unipotent_matrix
export partial_shift_graph
export partial_shift_graph_vertices
export rothe_matrix
export uniform_hypergraph

function lex_min_basis(A::MatElem)
  rref!(A)
  col_indices = Int[]
  row_index = 1
  row_bound = nrows(A)
  for col_index in 1:ncols(A)
    if !is_zero(A[row_index, col_index])
      push!(col_indices, col_index)
      row_index += 1
    end
    row_index > row_bound && break
  end
  return col_indices
end

# greedily find lex min basis for column span
function lex_min_basis(M::MatElem{<:MPolyRingElem})
  basis_indices = Int[]
  row_pivot_pos = 0
  for j in 1:ncols(M)
    # divide col by gcd of its entries
    c = content(M[:, [j]])
    if !isone(c)
      M[:, j] = divexact(M[:, [j]], c)
    end

    best_i = 0
    min_terms = typemax(Int)
    # find best entry in the col for the pivot
    # units are best, otherwise take entry with smallest number of terms
    for i in row_pivot_pos + 1:nrows(M)
      is_zero_entry(M, i, j) && continue
      if is_unit(M[i, j])
        # might not be best unit to take
        best_i = i
        break
      end
      # find row entry with smallest number of terms
      if length(M[i, j]) < min_terms
        min_terms = length(M[i, j])
        best_i = i
      end
    end

    best_i == 0 && break # column is lin comb of prev cols
    row_pivot_pos += 1
    push!(basis_indices, j)
    row_pivot_pos > nrows(M) && break # found all pivots
    M = swap_rows!(M, row_pivot_pos, best_i)
    j == ncols(M) && break
    for ii in 1:row_pivot_pos - 1
      jj = basis_indices[ii]
      _, a, b = gcd_with_cofactors(M[ii, j + 1], M[ii, jj])
      M[:, j + 1] = b * M[:, [j + 1]] - a * M[:, [jj]]
      M[:, j + 1] = divexact(M[:, [j + 1]], content(M[:, [j + 1]])) # divide col by gcd
      println(M[:, j + 1])
    end
  end
  return basis_indices
end

@doc raw"""
   symmetric_shift(F::Field, K::SimplicialComplex)

Returns the symmetric shift of `K`

see survey on algebraic shifting Gil Kalai (properly cite at some point)
"""
function symmetric_shift(F::Field, K::SimplicialComplex;
                         change_of_basis::T = nothing) where T <: Union{Nothing, MatElem}

  n = n_vertices(K)
  is_generic = false
  if isnothing(change_of_basis)
    is_generic = true
    # generic change of basis
    Fy, y = polynomial_ring(F, :y => (1:n, 1:n))
    change_of_basis = matrix(Fy, hcat(y))
    matrix_base = Fy
  else
    matrix_base = base_ring(change_of_basis)
    if matrix_base isa MPolyRing
      is_generic = true
    end
  end

  Rx, x = polynomial_ring(matrix_base, n)  
  # the generators of the stanley reisner ideal are combinations of [x_1, ..., x_n]
  R_K, _ = stanley_reisner_ring(Rx, K)

  # the computation is a over the field of fractions Fyx
  # we use a different ring to generate monomial_basis, coefficients need to be a field,
  # but we want to avoid using fraction field of Ry during row reduction 
  mb_ring, z = graded_polynomial_ring(F, n)


  input_faces = apply_on_faces(K) do (dim_face, faces)
    r = dim_face + 1

    mb = reverse(monomial_basis(mb_ring, r))
    A = Vector{elem_type(matrix_base)}[]
    mb_exponents = first.(collect.(exponents.(mb))) # gets monomial exponents

    for b in mb
      # need to compare with some alternatives
      transformed_monomial = evaluate(b, change_of_basis * gens(R_K))

      # we need to iterate through exponents here since the functions terms, coefficients or exponents
      # do not return 0 terms and we need to make sure the generic col aligns with the others
      # this is needed for the case when the field has finite characteristic
      # we use the lift because currently there is no function for get the coeff of
      # a MPolyQuoRingElem, which means we also need to check if it's zero before adding the coefficient
      col = [
        !is_zero(R_K(monomial(Rx, e))) ? coeff(lift(transformed_monomial), e) : matrix_base(0)
        for e in mb_exponents]

      push!(A, col)
    end
    
    C = matrix(matrix_base, reduce(hcat, A))

    if is_generic
      Oscar.ModStdQt.ref_ff_rc!(C)
    else
      rref!(C)
    end
    
    smallest_basis_el = z[r]^r
    smallest_index = findfirst(a -> a == smallest_basis_el, mb)
    col_indices = filter(x -> x >= smallest_index, independent_columns(C))
    monomial_exponents = first.(exponents.(mb[col_indices]))

    shifted_sets = Vector{Int}[]
    for me in monomial_exponents
      shifted_set = Int[]
      index_count = 1
      for (i, e) in enumerate(me)
        for j in 1:e
          push!(shifted_set, i - (r - index_count))
          index_count += 1
        end
      end
      push!(shifted_sets, shifted_set)
    end
    return shifted_sets
  end

  return simplicial_complex(input_faces)
end
