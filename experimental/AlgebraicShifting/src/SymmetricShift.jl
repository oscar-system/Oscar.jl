###############################################################################
# Symmetric shift
###############################################################################

"""
    symmetric_shift(F::Field, K::SimplicialComplex)

Returns the symmetric shift of `K`

TODO: add more docs and and expose for users
"""
function symmetric_shift(F::Field, K::SimplicialComplex, p::PermGroupElem)
  n = n_vertices(K)
  Fy, y = polynomial_ring(F, :y => (1:n, 1:n); cached=false)
  change_of_basis = rothe_matrix(Fy, p) * permutation_matrix(ZZ, p) * generic_unipotent_matrix(Fy)
  
  Rx, x = polynomial_ring(Fy, n)
  # the generators of the stanley reisner ideal are combinations of [x_1, ..., x_n]
  R_K, _ = stanley_reisner_ring(Rx, K)

  # the computation is a over the field of fractions Fyx
  # we use a different ring to generate monomial_basis, coefficients need to be a field,
  # but we want to avoid using fraction field of Ry during row reduction
  mb_ring, z = graded_polynomial_ring(F, n; cached=false)
  input_faces = Vector{Int}[]
  for r in 2:dim(K) + 1
    mb = reverse(monomial_basis(mb_ring, r))
    A = Vector{elem_type(Fy)}[]
    mb_exponents = first.(collect.(exponents.(mb))) # gets monomial exponents

    for b in mb
      # need to compare with some alternatives
      transformed_monomial = evaluate(b, change_of_basis * gens(R_K))

      # we need to iterate through exponents here since the functions terms, coefficients or exponents
      # do not return 0 terms and we need to make sure the generic col aligns with the others
      # this is needed for the case when the field has finite characteristic
      # we use the lift because currently there is no function to get the coeff of
      # a MPolyQuoRingElem, which means we also need to check if it's zero before adding the coefficient
      col = [
        !is_zero(R_K(monomial(Rx, e))) ? coeff(lift(transformed_monomial), e) : Fy(0)
        for e in mb_exponents]
      push!(A, col)
    end
    C = matrix(Fy, reduce(hcat, A))
    Oscar.ModStdQt.ref_ff_rc!(C)
    smallest_basis_el = z[r]^r
    smallest_index = findfirst(a -> a == smallest_basis_el, mb)
    col_indices = filter(x -> x >= smallest_index, independent_columns(C))
    monomial_exponents = first.(exponents.(mb[col_indices]))

    # adjustment to convert monomial to simplex
    for me in monomial_exponents
      shifted_set = Int[]
      index_count = 1
      for (i, e) in enumerate(me)
        for j in 1:e
          push!(shifted_set, i - (r - index_count))
          index_count += 1
        end
      end
      push!(input_faces, shifted_set)
    end
  end

  return simplicial_complex(input_faces)
end

function symmetric_shift(F::Field, K::SimplicialComplex, w::WeylGroupElem)
  iso = isomorphism(PermGroup, parent(w))
  return symmetric_shift(F, K, iso(w))
end
