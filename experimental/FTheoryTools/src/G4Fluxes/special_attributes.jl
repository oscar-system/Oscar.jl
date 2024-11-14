@doc raw"""
    ambient_space_models_of_g4_fluxes(m::AbstractFTheoryModel; check::Bool = true)::Vector{CohomologyClass}

Given an F-theory model $m$ defined as hypersurface in a simplicial and
complete toric base, we this method first computes a basis of
$H^(2,2)(X, \mathbb{Q})$ (by use of the method `basis_of_h22` below) and then filters
out "some" basis elements whose restriction to the hypersurface in question
is trivial. The exact meaning of "some" is explained above this method.

Note that it can be computationally very demanding to check if a toric variety
$X$ is complete (and simplicial). The optional argument `check` can be set
to `false` to skip these tests.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> Kbar = anticanonical_divisor_class(B3)
Divisor class on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w"=>Kbar))
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> g4_amb_list = ambient_space_models_of_g4_fluxes(t)
2-element Vector{CohomologyClass}:
 Cohomology class on a normal toric variety given by z^2
 Cohomology class on a normal toric variety given by y^2

julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 8))
Hypersurface model over a concrete base

julia> g4_amb_list = ambient_space_models_of_g4_fluxes(qsm_model, check = false);

julia> length(g4_amb_list) == 172
true
```
"""
function ambient_space_models_of_g4_fluxes(m::AbstractFTheoryModel; check::Bool = true)::Vector{CohomologyClass}

  # Entry check
  @req base_space(m) isa NormalToricVariety "Base space must be a toric variety for computation of ambient space G4 candidates"
  if has_attribute(m, :ambient_space_models_of_g4_fluxes)
    return get_attribute(m, :ambient_space_models_of_g4_fluxes)
  end

  # Execute entry tests in computation of basis_of_h22. If any of these fail, no need to proceed. Hence, do this first.
  filtered_h22_basis = basis_of_h22(ambient_space(m), check = check)

  # Prepare data of the toric ambient space
  gS = gens(cox_ring(ambient_space(m)))
  mnf = Oscar._minimal_nonfaces(ambient_space(m))
  sr_ideal_pos = [Vector{Int}(Polymake.row(mnf, i)) for i in 1:Polymake.nrows(mnf)]

  # Filter out basis elements
  for a in length(filtered_h22_basis):-1:1

    # Find non-zero exponent positions in the polynomial filtered_h22_basis[a]
    exp_list = collect(exponents(polynomial(filtered_h22_basis[a]).f))[1]
    vanishing_vars_pos = findall(!=(0), exp_list)
    
    # Simplify the hypersurface polynomial by setting relevant variables to zero
    new_pt = divrem(hypersurface_equation(m), gS[vanishing_vars_pos[1]])[2]
    if length(vanishing_vars_pos) == 2
      new_pt = divrem(new_pt, gS[vanishing_vars_pos[2]])[2]
    end

    # If all coefficient of `new_pt` sum to zero, keep this generator.
    if sum(coefficients(new_pt)) == 0
      continue
    end
    
    # Determine remaining variables, after scaling "away" others.
    remaining_vars_list = Set(1:length(gS))
    for my_exps in sr_ideal_pos
      len_my_exps = length(my_exps)
      inter_len = count(idx -> idx in vanishing_vars_pos, my_exps)
      if (len_my_exps == 2 && inter_len == 1) || (len_my_exps == 3 && inter_len == 2)
        delete!(remaining_vars_list, my_exps[findfirst(idx -> !(idx in vanishing_vars_pos), my_exps)])
      end
    end
    remaining_vars_list = collect(remaining_vars_list)

    # If one monomial of `new_pt` has unset positions, then keep this generator.
    delete_it = true
    for exps in exponents(new_pt)
      if any(x -> x != 0, exps[remaining_vars_list])
        delete_it = false
        break
      end
    end
    if delete_it
      deleteat!(filtered_h22_basis, a)
    end

  end

  set_attribute!(m, :ambient_space_models_of_g4_fluxes, filtered_h22_basis)
  return filtered_h22_basis

end


@doc raw"""
    well_quantized_ambient_space_models_of_g4_fluxes(m::AbstractFTheoryModel; check::Bool = true)::Tuple{QQMatrix, QQMatrix}

Given an F-theory model $m$ defined as hypersurface in a simplicial and
complete toric base, this method computes a basis of all well-quantized
ambient space $G_4$-fluxes. The result of this operation is a tuple of two
matrices. The columns of the first matrix specify those (rational) combinations
of ambient space $G_4$-fluxes, of which one may only take $\mathbb{Z}$-linear
combinations without violating flux quantization. The columns of the second
matrix specify those (rational) combinations of ambient space $G_4$-fluxes, for
which any rational linear combination satisfies the flux quantization condition.

It can be computationally very demanding to check if a toric variety
$X$ is complete (and simplicial). The optional argument `check` can be set
to `false` to skip these tests.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> Kbar = anticanonical_divisor_class(B3)
Divisor class on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w"=>Kbar))
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> ambient_space_models_of_g4_fluxes(t, check = false);

julia> res = well_quantized_ambient_space_models_of_g4_fluxes(t, check = false);

julia> res[1]
[1//4   -3//16]
[   0   1//144]

julia> res[2]
2 by 0 empty matrix
```
"""
function well_quantized_ambient_space_models_of_g4_fluxes(m::AbstractFTheoryModel; check::Bool = true)::Tuple{QQMatrix, QQMatrix}

  # (1) Entry checks
  @req base_space(m) isa NormalToricVariety "Computation of well-quantized G4-fluxes only supported for toric base and ambient spaces"
  @req dim(ambient_space(m)) == 5 "Computation of well-quantized G4-fluxes only supported for 5-dimensional toric ambient spaces"
  if check
    @req is_complete(ambient_space(m)) "Computation of well-quantized G4-fluxes only supported for complete toric ambient spaces"
    @req is_simplicial(ambient_space(m)) "Computation of well-quantized G4-fluxes only supported for simplicial toric ambient space"
  end
  if has_attribute(m, :well_quantized_ambient_space_models_of_g4_fluxes)
    return get_attribute(m, :well_quantized_ambient_space_models_of_g4_fluxes)
  end


  # (2) Compute data, that is frequently used by the sophisticated intersection product below
  S = cox_ring(ambient_space(m))
  gS = gens(cox_ring(ambient_space(m)))
  linear_relations = matrix(QQ, matrix(ZZ, rays(ambient_space(m))))
  scalings = [c.coeff for c in S.d]
  mnf = Oscar._minimal_nonfaces(ambient_space(m))
  sr_ideal_pos = [Vector{Int}(Polymake.row(mnf, i)) for i in 1:Polymake.nrows(mnf)]
  data = (
    S = S,
    gS = gS,
    linear_relations = linear_relations,
    scalings = scalings,
    sr_ideal_pos = sr_ideal_pos
  )


  # (3) Are intersection numbers known?
  inter_dict = Dict{NTuple{4, Int64}, QQFieldElem}()
  s_inter_dict = Dict{String, QQFieldElem}()
  if has_attribute(m, :inter_dict)
    inter_dict = get_attribute(m, :inter_dict)
  end
  if has_attribute(m, :s_inter_dict)
    s_inter_dict = get_attribute(m, :s_inter_dict)
  end


  # (4) Obtain critical information - this may take significant time!
  ambient_space_flux_candidates_basis = ambient_space_models_of_g4_fluxes(m, check = check)
  list_of_divisor_pairs_to_be_considered = Oscar._ambient_space_divisors_to_be_considered(m)


  # (5) Work out the relevant intersection numbers and organize them in a constraint_matrix.
  constraint_matrix = Vector{Vector{QQFieldElem}}()
  if arxiv_doi(m) == "10.48550/arXiv.1511.03209"

    # Use special intersection theory for special F-theory model. This technology could be extended beyond this one use-case in the future.
    for i in 1:length(ambient_space_flux_candidates_basis)
      condition = Vector{ZZRingElem}()
      idxs = findall(x -> x != 0, collect(exponents(polynomial(ambient_space_flux_candidates_basis[i]).f))[1])
      if length(idxs) == 1
        idxs = [idxs[1], idxs[1]]
      end

      for j in 1:length(list_of_divisor_pairs_to_be_considered)

        my_tuple = Tuple(sort([idxs..., list_of_divisor_pairs_to_be_considered[j]...]))
        push!(condition, sophisticated_intersection_product(ambient_space(m), my_tuple, hypersurface_equation(m), inter_dict, s_inter_dict, data))

      end
      push!(constraint_matrix, condition)
    end

  else
  
    # Cover all other case with generic, but potentially painfully slow methodology.
    tds = torusinvariant_prime_divisors(ambient_space(m))
    cds = [cohomology_class(tds[i]) for i in 1:length(tds)]
    pt_class = cohomology_class(anticanonical_divisor_class(ambient_space(m)))
    for i in 1:length(ambient_space_flux_candidates_basis)

      condition = Vector{QQFieldElem}()

      idxs = findall(x -> x != 0, collect(exponents(polynomial(ambient_space_flux_candidates_basis[i]).f))[1])
      if length(idxs) == 1
        idxs = [idxs[1], idxs[1]]
      end

      for j in 1:length(list_of_divisor_pairs_to_be_considered)

        my_tuple = Tuple(sort([idxs..., list_of_divisor_pairs_to_be_considered[j]...]))
        if !haskey(inter_dict, my_tuple)
          class = ambient_space_flux_candidates_basis[i] * cds[list_of_divisor_pairs_to_be_considered[j][1]] * cds[list_of_divisor_pairs_to_be_considered[j][2]] * pt_class
          inter_dict[my_tuple] = integrate(class)
        end
        push!(condition, inter_dict[my_tuple])

      end
      push!(constraint_matrix, condition)
    end

  end
  
  # (6) Convert the intersection matrix to a ZZ matrix. If necessary, multiply it by a suitable integer.
  # (6) Then compute its Smith normal form.
  denom = lcm(unique(vcat([denominator.(k) for k in constraint_matrix]...)))
  if denom != 1
    constraint_matrix = denom * constraint_matrix
  end
  C = transpose(matrix(ZZ, constraint_matrix))
  S, T, U = snf_with_transform(C)


  # (7) Recall that we are seeking constraint_matrix * q = n, where q is a vector of rational numbers and n a vector of natural numbers.
  # (7) Above, we multiplied the constraint matrix with the denominator if necessary. Thereby, this equation is equivalent to C * q = d * n.
  # (7) Since T is an invertible matrix with integer entries, this system is equivalent to T * C * q = T * (d * n).
  # (7) This in turn is equivalent to (T * C * U) * (U^-1 * q) = d * (T*n).
  # (7) In other words, we have S * (U^(-1) * q) = d * (T*n).
  # (7) Note that T*n are just other integers, which we may use to parametrize the space of solutions. Let us thus write:
  # (7) So we have S * (U^-1) * q) = d * N, with N the newly chosen parametrization. Let tilde_q = U^(-1) * q.
  r = rank(S)
  @req all(k -> !is_zero(S[k,k]), 1:r) "Inconsistency in Smith normal form detected. Please inform the authors."
  @req all(k -> is_zero(S[k,k]), r+1:min(nrows(S), ncols(S))) "Inconsistency in Smith normal form  detected. Please inform the authors."
  # (7) S is diagonal, and has non-zero entries at diagonal position 1 to r: S = (l1, ..., lr, 0, ..., 0). Therefore, every q_tilde
  # (7) which solves the above is of the form q_tilde = (1/l1 * d * N1, ..., 1/lr * d * N_r, tilde_q_(r+1), ..., ), where 
  # (7) tilde_q_(r+1) etc. are unconstrained. We encode this solution in the following matrix.
  S_prime = zero_matrix(QQ, ncols(S), ncols(S))
  for k in 1:min(nrows(S), ncols(S))
    if k <= r
      S_prime[k,k] = denom//S[k,k]
    else
      S_prime[k,k] = 1
    end
  end
  # (7) To extract the solutions q from these solutions tilde_q, we multiply with the matrix U from the left.
  solution_matrix = U * S_prime


  # (8) Overall, we are allowed to take any Z-linear combinations of the first r columns of solution_matrix together with any
  # (8) rational combination of its remaining columns. Those are exactly the G4-fluxes which pass the elementary quantization tests.
  # (8) Note that the set of vectors for which we can allow any rational combination is isomorphic to the kernel of constraint_matrix.
  # (8) Indeed, only this ensures that upon multiplication with any rational number, the result remains an integer.
  res = (solution_matrix[:,1:r], solution_matrix[:,r+1:ncols(solution_matrix)])


  # (9) Remember computed data
  set_attribute!(m, :well_quantized_ambient_space_models_of_g4_fluxes, res)
  set_attribute!(m, :inter_dict, inter_dict)
  set_attribute!(m, :s_inter_dict, s_inter_dict)


  # (10) Finally, return the result
  return res

end
