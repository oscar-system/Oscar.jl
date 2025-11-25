
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
    g = matrix(Fy, reduce(hcat, A))
    append!(input_faces, symmetric_shift(K, g))
  end

  return simplicial_complex(input_faces)
end

function symmetric_shift(F::Field, K::SimplicialComplex, w::WeylGroupElem)
  iso = isomorphism(PermGroup, parent(w))
  return symmetric_shift(F, K, iso(w))
end

function symmetric_shift(K::UniformHypergraph, g::MatElem)
  faces = Vector{Int}[]
  Oscar.ModStdQt.ref_ff_rc!(g)
  smallest_basis_el = z[r]^r
  smallest_index = findfirst(a -> a == smallest_basis_el, mb)
  col_indices = filter(x -> x >= smallest_index, pivot_columns(g))
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
    push!(faces, shifted_set)
  end
  return faces
end

function symmetric_shift_lv(F::Field, K::ComplexOrHypergraph, p::PermGroupElem; n_samples=100, timed=false, kw...)
  # this might need to be changed based on the characteristic
  # we expect that the larger the characteristic the smaller the sample needs to be
  # setting to 100 now for good measure
  # Compute n_samples many shifts by radom matrices, and take the lexicographically minimal one, together with its first index of occurrence.
  (shift, i), stats... = @timed efindmin((exterior_shift(K, random_rep_bruhat_cell(F, p); kw...) for i in 1:n_samples); lt=isless_lex)
  # Check if `shift` is the generic exterior shift of K
  prime_field = characteristic(F) == 0 ? QQ : fpField(UInt(characteristic(F)))
  n = n_vertices(K)
  is_correct_shift, stats2... = @timed (p != perm(reverse(1:n)) || is_shifted(shift)) && check_shifted(prime_field, K, shift, p; kw...)
  
  if timed
    if is_correct_shift
      return shift, (i, stats.time, stats.bytes, stats2.time, stats2.bytes)
    else
      return nothing, (">$n_samples", stats.time, stats.bytes, stats2.time, stats2.bytes)
    end
  else
    return is_correct_shift ? shift : nothing
  end
end
