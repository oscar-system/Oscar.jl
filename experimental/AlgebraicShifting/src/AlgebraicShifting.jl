function independent_columns(A::MatElem)
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

# applies function on faces of equal dimension
function apply_on_faces(f::Function, K::SimplicialComplex)
  dim_K = dim(K)
  K_facets = facets(K)
  facet_dict = Dict(
    (k => Set(filter(f -> length(f) == k + 1, K_facets)) for k in 1:dim_K )
  )

  # look into hasse diagram from polymake
  
  for facet in K_facets
    dim_facet = length(facet) - 1
    
    for subset_size in 2:dim_facet
      for subset in subsets(facet, subset_size)
        push!(facet_dict[subset_size - 1], subset)
      end
    end
  end

  input_faces = Vector{Int}[]
  for (dim_face, faces) in facet_dict
    input_faces = vcat(input_faces, f((dim_face, collect.(faces))))
  end
  return input_faces
end

function shift_basis_ext(F::Field, K::SimplicialComplex)
  n = n_vertices(K)
  shifted_K = Oscar.delta_ext(F, K)
  lambda, _ = exterior_algebra(F, n)
  x = gens(lambda)
  
  ker_gens = [
    reduce(*, lift.(
      x[collect(f)]
    )) for f in minimal_nonfaces(K)
      ]

  lift_ring = base_ring(lambda)
  J = two_sided_ideal(lift_ring, (ker_gens))
  A, _ = quo(lift_ring, J + modulus(lambda))
  x = gens(A)
  basis = [reduce(*, x[collect(s)]) for s in facets(shifted_K)]

  return A, basis
end

@doc raw"""
     exterior_shift(F::Field, K::SimplicialComplex)

Returns the exterior shift of `K`
"""
function exterior_shift(F::Field, K::SimplicialComplex;
                        change_of_basis::T = nothing) where T <: Union{Nothing, MatElem}
  # the exterior shifting works in a different algebra that lends
  # itself to an easier implementation 
  n = n_vertices(K)
  is_generic = false
  if isnothing(change_of_basis)
    is_generic = true
    # generic change of basis
    Fx, x = polynomial_ring(F, :x => (1:n, 1:n))
    change_of_basis = matrix(Fx, hcat(x))
    matrix_base = Fx
  else
    matrix_base = base_ring(change_of_basis)
    if matrix_base isa MPolyRing
      is_generic = true
    end
  end

  # might not be needed
  cmp(S1, S2) = min(symdiff(S1, S2)...) in S1

  input_faces = apply_on_faces(K) do (dim_face, faces)
    # is needed here since otherwise the sets are not
    # in lex order, seubsets is being called from Hecke,
    # maybe we should look into this and see if polymake's
    # is preffered?
    nCk = sort(subsets(n, dim_face + 1))
    entry_type = elem_type(matrix_base)

    sub_compound_matrix = Vector{entry_type}[]

    for col_subset in nCk
      row_minors = entry_type[]
      for row_subset in faces
        ri = collect(row_subset)
        ci = collect(col_subset)
        push!(row_minors, det(change_of_basis[ri, ci]))
      end
      push!(sub_compound_matrix, row_minors)
    end
    A = matrix(reduce(hcat, sub_compound_matrix))

    if is_generic
      Oscar.ModStdQt.ref_ff_rc!(A)
    else
      rref!(A)
    end
    delta_k = independent_columns(A)
    return nCk[delta_k]
  end

  return simplicial_complex(input_faces)
end

@doc raw"""
   symmetric_shift(F::Field, K::SimplicialComplex)

Returns the symmetric shift of `K`

see survey on algebraic shifting Gil Kalai (properly cite at some point)
"""
function symmetric_shift(F::Field, K::SimplicialComplex)
  n = n_vertices(K)
  Fy, y = polynomial_ring(F, :y => (1:n, 1:n))
  Fyx, x = polynomial_ring(Fy, n)
  # the computation is a over the field of fractions Fyx
  # we use a different ring to generate monomial_basis, coefficients need to be a field,
  # but we want to avoid using fraction field of Ry during row reduction 
  mb_ring, z = graded_polynomial_ring(F, n)
  o = monomial_ordering(mb_ring, :lex)

  # the generators of the stanley reisner ideal are combinations of [x_1, ..., x_n]
  R_K, _ = stanley_reisner_ring(Fyx, K)

  input_faces = apply_on_faces(K) do (dim_face, faces)
    r = dim_face + 1

    mb = reverse(monomial_basis(mb_ring, r))
    Y = matrix(Fy, hcat(y))
    A = Vector{MPolyRingElem}[]

    for b in mb
      # need to compare with some alternatives
      transformed_monomial = evaluate(b, Y * gens(R_K))

      # double check this line
      # is stanley reisner ring over Z needed to then reduce later
      # stanley reisner over F is cheaper
      # this part will need to be adjusted for finite fields
      # there should be no lifting
      non_zero_terms = filter(x -> !is_zero(R_K(x)), collect(terms(lift(transformed_monomial))))

      
      # coefficients is a function for polynomials but since we are looking at the terms
      # i.e. we are seperating by monomials already, we know the first coefficient is just
      # the coefficient of the monomial
      generic_col = first.(coefficients.(non_zero_terms))
      push!(A, generic_col)
    end
    C = matrix(Fy, reduce(hcat, A))
    Oscar.ModStdQt.ref_ff_rc!(C)
    
    smallest_basis_el = z[r]^r
    smallest_index = findfirst(a -> a == smallest_basis_el, mb)
    col_indices = filter(x -> x >= smallest_index, independent_columns(C))
    generic_monomial_exponents = first.(exponents.(mb[col_indices]))

    shifted_sets = Vector{Int}[]
    for me in generic_monomial_exponents
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
