@doc raw"""
    special_flux_family(m::AbstractFTheoryModel; not_breaking::Bool = false, check::Bool = true)

Compute a family of G4-fluxes with specified properties for a given F-theory model `m`,
defined as a hypersurface in a simplicial and complete toric ambient space.

### Description
This method models the G4-flux family using the restriction of cohomology classes on the
toric ambient space to the hypersurface. In the toric ambient space, these classes are vertical,
meaning they are of the form $a \wedge b$, where $a, b \in H^(1,1)(X_\Sigma)$, with $X_\Sigma$
denoting the toric ambient space.

The resulting family is subjected to consistency conditions to ensure it satisfies elementary
quantization requirements and transversality conditions. By default, the method returns ambient
space G4-flux candidates that meet these conditions.

### Optional Arguments
- `not_breaking = true`: Ensures the flux family preserves the non-abelian gauge group.
- `check = false`: Skips computational checks for whether the toric ambient space $X_\Sigma$ is complete and simplicial, which can be resource-intensive.

### Notes
!!! warning
    This method assumes that $c_2( \widehat{Y}_4)$ is even. No checks or errors are implemented for this condition, so use cautiously.

This assumption relates to the quantization condition, which requires verifying if the twist of a
given $G_4$-flux by $1/2 \cdot c_2( \widehat{Y}_4)$ is even. For many F-theory models, such as
compactifications on smooth Calabi-Yau 4-folds with globally defined Weierstrass models ([CS12](@cite)),
$c_2( \widehat{Y}_4)$ is known to be even. This also applies to all F-theory QSMs ([CHLLT19](@cite)).

### Computational Details
The method internally identifies two matrices related to the family of fluxes:
1. `matrix_integral`: Specifies rational combinations of ambient space G4-flux candidates that can only form $\mathbb{Z}$-linear combinations without violating elementary flux quantization conditions.
2. `matrix_rational`: Specifies rational combinations of flux candidates for which any rational linear combination satisfies the elementary flux quantization conditions.

These matrices are accessible for further analysis.

### Examples

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> fg = special_flux_family(qsm_model, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: broken

julia> g4_tester = random_flux_instance(fg, check = false)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not executed

julia> g4_tester_double = g4_flux(qsm_model, cohomology_class(g4_tester), check = false);

julia> is_well_quantized(g4_tester_double)
true

julia> passes_transversality_checks(g4_tester_double)
true

julia> c = [60, 51, 90, 0, 24, 51, -24, 45, 30, 0, -48, 90, -57, 60, 30, 15, 120, 0, -60, 0, -720, -420, -270, -60, -2190];

julia> qsm_g4_flux = flux_instance(fg, transpose(matrix(ZZ, [c])), zero_matrix(QQ, 0, 1), check = false)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not executed

julia> passes_transversality_checks(qsm_g4_flux)
true

julia> is_well_quantized(qsm_g4_flux)
true
```

Finally, we demonstrate the computation of the well-quantized $G_4$-fluxes which pass the transversality checks,
and which in addition do not break the non-abelian gauge group.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> fg = special_flux_family(qsm_model, not_breaking = true, check = false)
A family of G4 fluxes:
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: not broken

julia> g4_tester = random_flux_instance(fg, check = false)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: not broken
  - Tadpole cancellation check: not executed

julia> g4_tester_double = g4_flux(qsm_model, cohomology_class(g4_tester), check = false);

julia> is_well_quantized(g4_tester_double)
true

julia> passes_transversality_checks(g4_tester_double)
true

julia> breaks_non_abelian_gauge_group(g4_tester_double)
false

julia> c = [3];

julia> qsm_g4_flux = flux_instance(fg, matrix(ZZ, [[3]]), zero_matrix(QQ, 0, 1), check = false)
G4-flux candidate
  - Elementary quantization checks: satisfied
  - Transversality checks: satisfied
  - Non-abelian gauge group: not broken
  - Tadpole cancellation check: not executed

julia> is_well_quantized(qsm_g4_flux)
true

julia> passes_transversality_checks(qsm_g4_flux)
true

julia> breaks_non_abelian_gauge_group(qsm_g4_flux)
false

julia> divs = torusinvariant_prime_divisors(ambient_space(qsm_model));

julia> e1 = cohomology_class(divs[15]);

julia> e2 = cohomology_class(divs[12]);

julia> e4 = cohomology_class(divs[14]);

julia> u = cohomology_class(divs[13]);

julia> v = cohomology_class(divs[10]);

julia> pb_Kbar = cohomology_class(sum([divs[k] for k in 1:9]));

julia> g4_class = (-3) // kbar3(qsm_model) * (5 * e1 * e4 + pb_Kbar * (-3 * e1 - 2 * e2 - 6 * e4 + pb_Kbar - 4 * u + v));

julia> qsm_g4_flux == g4_flux(qsm_model, g4_class)
true
```
"""
function special_flux_family(m::AbstractFTheoryModel; not_breaking::Bool = false, check::Bool = true)
  if !not_breaking
    return well_quantized_and_transversal_ambient_space_models_of_g4_fluxes(m, check = check)
  else
    return well_quantized_and_transversal_and_no_non_abelian_gauge_group_breaking_ambient_space_models_of_g4_fluxes(m, check = check)
  end
end


@attr FamilyOfG4Fluxes function well_quantized_ambient_space_models_of_g4_fluxes(m::AbstractFTheoryModel; check::Bool = true)

  # (0) Has this result been computed before?
  if has_attribute(m, :matrix_integral_quant) && has_attribute(m, :matrix_rational_quant)
    fgs = family_of_g4_fluxes(m, matrix_integral_quant(m, check = check), matrix_rational_quant(m, check = check))
    set_attribute!(fgs, :is_well_quantized, true)
    set_attribute!(fgs, :passes_transversality_checks, false)
    set_attribute!(fgs, :breaks_non_abelian_gauge_group, true)
    return fgs
  end


  # (1) Entry checks
  @req base_space(m) isa NormalToricVariety "Computation of well-quantized G4-fluxes only supported for toric base and ambient spaces"
  @req dim(ambient_space(m)) == 5 "Computation of well-quantized G4-fluxes only supported for 5-dimensional toric ambient spaces"
  if check
    @req is_complete(ambient_space(m)) "Computation of well-quantized G4-fluxes only supported for complete toric ambient spaces"
    @req is_simplicial(ambient_space(m)) "Computation of well-quantized G4-fluxes only supported for simplicial toric ambient space"
  end


  # (2) Compute data, that is frequently used by the sophisticated intersection product below
  S = cox_ring(ambient_space(m))
  gS = gens(cox_ring(ambient_space(m)))
  linear_relations = matrix(QQ, rays(ambient_space(m)))
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
  # TODO: If available and necessary, convert inter_dict.
  # TODO: This is necessary, because serializing and loading turns NTuple{4, Int64} into Tuple (as of March 5, 2025).
  # TODO: Once serialization has caught up, this conversion will no longer be needed.
  if has_attribute(m, :inter_dict) && typeof(get_attribute(m, :inter_dict)) != Dict{NTuple{4, Int64}, ZZRingElem}
    original_dict = get_attribute(m, :inter_dict)
    new_dict = Dict{NTuple{4, Int64}, ZZRingElem}()
    for (key, value) in original_dict
      new_key = NTuple{4, Int64}(key)
      new_dict[new_key] = value
    end
    set_attribute!(model, :inter_dict, new_dict)
  end
  inter_dict = get_attribute!(m, :inter_dict) do
    Dict{NTuple{4, Int64}, ZZRingElem}()
  end::Dict{NTuple{4, Int64}, ZZRingElem}
  s_inter_dict = get_attribute!(m, :s_inter_dict) do
    Dict{String, ZZRingElem}()
  end::Dict{String, ZZRingElem}


  # (4) Obtain critical information - this may take significant time!
  ambient_space_flux_candidates_basis = _ambient_space_models_of_g4_fluxes(m, check = check)
  ambient_space_flux_candidates_basis_indices = get_attribute(m, :ambient_space_models_of_g4_fluxes_indices)::Vector{Tuple{Int64, Int64}}
  list_of_divisor_pairs_to_be_considered = Oscar._ambient_space_divisor_pairs_to_be_considered(m)


  # (5) Work out the relevant intersection numbers and organize them in a constraint_matrix.
  constraint_matrix = Vector{Vector{QQFieldElem}}()
  # I have prepared some functionality below, regarding the case that this matrix should have rational entries.
  # However, I expect that this will not happen as long as the hypersurface in question is smooth.
  if arxiv_doi(m) == "10.48550/arXiv.1511.03209"

    # Use special intersection theory for special F-theory model. This technology could be extended beyond this one use-case in the future.
    for i in 1:length(ambient_space_flux_candidates_basis)
      condition = Vector{ZZRingElem}()
      for j in 1:length(list_of_divisor_pairs_to_be_considered)
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., list_of_divisor_pairs_to_be_considered[j]...]))
        push!(condition, sophisticated_intersection_product(ambient_space(m), my_tuple, hypersurface_equation(m), inter_dict, s_inter_dict, data))
      end
      push!(constraint_matrix, condition)
    end

  else
  
    # Cover all other case with generic, but potentially painfully slow methodology.
    tds = torusinvariant_prime_divisors(ambient_space(m))
    cds = [cohomology_class(td) for td in tds]
    pt_class = cohomology_class(anticanonical_divisor_class(ambient_space(m)))
    for i in 1:length(ambient_space_flux_candidates_basis)
      condition = Vector{ZZRingElem}()
      for j in 1:length(list_of_divisor_pairs_to_be_considered)
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., list_of_divisor_pairs_to_be_considered[j]...]))
        push!(condition, get!(inter_dict, my_tuple) do
          class = ambient_space_flux_candidates_basis[i] * cds[list_of_divisor_pairs_to_be_considered[j][1]] * cds[list_of_divisor_pairs_to_be_considered[j][2]] * pt_class
          return ZZ(integrate(class))
        end)
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
  fgs = family_of_g4_fluxes(m, res[1], res[2])
  set_attribute!(m, :matrix_integral_quant, res[1])
  set_attribute!(m, :matrix_rational_quant, res[2])
  set_attribute!(fgs, :is_well_quantized, true)
  set_attribute!(fgs, :passes_transversality_checks, false)
  set_attribute!(fgs, :breaks_non_abelian_gauge_group, true)
  set_attribute!(m, :inter_dict, inter_dict)
  set_attribute!(m, :s_inter_dict, s_inter_dict)


  # (10) Finally, return the result
  return fgs
end


@attr FamilyOfG4Fluxes function well_quantized_and_transversal_ambient_space_models_of_g4_fluxes(m::AbstractFTheoryModel; check::Bool = true)

  # (0) Has this result been computed before?
  if has_attribute(m, :matrix_integral_quant_transverse) && has_attribute(m, :matrix_rational_quant_transverse)
    fgs = family_of_g4_fluxes(m, matrix_integral_quant_transverse(m, check = check), matrix_rational_quant_transverse(m, check = check))
    set_attribute!(fgs, :is_well_quantized, true)
    set_attribute!(fgs, :passes_transversality_checks, true)
    set_attribute!(fgs, :breaks_non_abelian_gauge_group, true)
    return fgs
  end
  
  
  # (1) Entry checks
  @req base_space(m) isa NormalToricVariety "Computation of well-quantized and transversal G4-fluxes only supported for toric base and ambient spaces"
  @req dim(ambient_space(m)) == 5 "Computation of well-quantized and transversal G4-fluxes only supported for 5-dimensional toric ambient spaces"
  if check
    @req is_complete(ambient_space(m)) "Computation of well-quantized and transversal G4-fluxes only supported for complete toric ambient spaces"
    @req is_simplicial(ambient_space(m)) "Computation of well-quantized and transversal G4-fluxes only supported for simplicial toric ambient space"
  end


  # (2) Compute data, that is frequently used by the sophisticated intersection product below
  S = cox_ring(ambient_space(m))
  gS = gens(cox_ring(ambient_space(m)))
  linear_relations = matrix(QQ, rays(ambient_space(m)))
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
  # TODO: If available and necessary, convert inter_dict.
  # TODO: This is necessary, because serializing and loading turns NTuple{4, Int64} into Tuple (as of March 5, 2025).
  # TODO: Once serialization has caught up, this conversion will no longer be needed.
  if has_attribute(m, :inter_dict) && typeof(get_attribute(m, :inter_dict)) != Dict{NTuple{4, Int64}, ZZRingElem}
    original_dict = get_attribute(m, :inter_dict)
    new_dict = Dict{NTuple{4, Int64}, ZZRingElem}()
    for (key, value) in original_dict
      new_key = NTuple{4, Int64}(key)
      new_dict[new_key] = value
    end
    set_attribute!(model, :inter_dict, new_dict)
  end
  inter_dict = get_attribute!(m, :inter_dict) do
    Dict{NTuple{4, Int64}, ZZRingElem}()
  end::Dict{NTuple{4, Int64}, ZZRingElem}
  s_inter_dict = get_attribute!(m, :s_inter_dict) do
    Dict{String, ZZRingElem}()
  end::Dict{String, ZZRingElem}


  # (4) Obtain critical information - this may take significant time!
  ambient_space_flux_candidates_basis = _ambient_space_models_of_g4_fluxes(m, check = check)
  list_of_base_divisor_pairs_to_be_considered = Oscar._ambient_space_base_divisor_pairs_to_be_considered(m)
  ambient_space_flux_candidates_basis_indices = get_attribute(m, :ambient_space_models_of_g4_fluxes_indices)::Vector{Tuple{Int64, Int64}}
  list_of_divisor_pairs_to_be_considered = Oscar._ambient_space_divisor_pairs_to_be_considered(m)


  # (5) Work out the relevant intersection numbers to tell if a flux is vertical
  vertical_constraint_matrix = Vector{Vector{QQFieldElem}}()
  # I have prepared some functionality below, regarding the case that this matrix should have rational entries.
  # However, I expect that this will not happen as long as the hypersurface in question is smooth.
  if arxiv_doi(m) == "10.48550/arXiv.1511.03209"

    # Use special intersection theory for special F-theory model. This technology could be extended beyond this one use-case in the future.
    for i in 1:length(ambient_space_flux_candidates_basis)

      condition = Vector{ZZRingElem}()

      # Compute against pairs of base divisors
      for j in 1:length(list_of_base_divisor_pairs_to_be_considered)
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., list_of_base_divisor_pairs_to_be_considered[j]...]))
        push!(condition, sophisticated_intersection_product(ambient_space(m), my_tuple, hypersurface_equation(m), inter_dict, s_inter_dict, data))
      end

      # Compute against zero section and base divisor
      pos_zero_section = findfirst(x -> x == "z", string.(gS))
      for j in 1:n_rays(base_space(m))
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., [j, pos_zero_section]...]))
        push!(condition, sophisticated_intersection_product(ambient_space(m), my_tuple, hypersurface_equation(m), inter_dict, s_inter_dict, data))
      end

      push!(vertical_constraint_matrix, condition)

    end

  else
  
    # Cover all other case with generic, but potentially painfully slow methodology.
    tds = torusinvariant_prime_divisors(ambient_space(m))
    cds = [cohomology_class(td) for td in tds]
    pt_class = cohomology_class(anticanonical_divisor_class(ambient_space(m)))
    for i in 1:length(ambient_space_flux_candidates_basis)

      condition = Vector{ZZRingElem}()

      # Compute against pairs of base divisors
      for j in 1:length(list_of_base_divisor_pairs_to_be_considered)
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., list_of_base_divisor_pairs_to_be_considered[j]...]))
        push!(condition, get!(inter_dict, my_tuple) do
          class = ambient_space_flux_candidates_basis[i] * cds[list_of_base_divisor_pairs_to_be_considered[j][1]] * cds[list_of_base_divisor_pairs_to_be_considered[j][2]] * pt_class
          return ZZ(integrate(class))
        end)
      end

      # Compute against zero section and base divisor
      zsc = zero_section_class(m)
      pos_zero_section = findfirst(x -> x == string(polynomial(zsc)), string.([polynomial(x) for x in cds]))
      @req pos_zero_section !== nothing && pos_zero_section >= 1 "Could not establish position of the zero section"
      for j in 1:n_rays(base_space(m))
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., [j, pos_zero_section]...]))
        push!(condition, get!(inter_dict, my_tuple) do
          class = ambient_space_flux_candidates_basis[i] * cds[j] * zsc * pt_class
          return ZZ(integrate(class))
        end)

      end

      push!(vertical_constraint_matrix, condition)

    end

  end


  # (6) Compute the vertical fluxes as the kernel of the vertical_constraint_matrix.
  # (6) To later tell if those fluxes are properly quantized, we want to parametrize them with integer coefficient only.
  denom = lcm(unique(vcat([denominator.(k) for k in vertical_constraint_matrix]...)))
  if denom != 1
    vertical_constraint_matrix = denom * vertical_constraint_matrix
  end
  C_vertical = transpose(matrix(ZZ, vertical_constraint_matrix))
  vertical_fluxes = nullspace(C_vertical)[2]


  # (7) Work out the relevant intersection numbers to tell if a flux is well quantized
  quant_constraint_matrix = Vector{Vector{QQFieldElem}}()
  # I have prepared some functionality below, regarding the case that this matrix should have rational entries.
  # However, I expect that this will not happen as long as the hypersurface in question is smooth.
  if arxiv_doi(m) == "10.48550/arXiv.1511.03209"

    # Use special intersection theory for special F-theory model. This technology could be extended beyond this one use-case in the future.
    for i in 1:length(ambient_space_flux_candidates_basis)
      condition = Vector{ZZRingElem}()
      for j in 1:length(list_of_divisor_pairs_to_be_considered)
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., list_of_divisor_pairs_to_be_considered[j]...]))
        push!(condition, sophisticated_intersection_product(ambient_space(m), my_tuple, hypersurface_equation(m), inter_dict, s_inter_dict, data))
      end
      push!(quant_constraint_matrix, condition)
    end

  else
  
    # Cover all other case with generic, but potentially painfully slow methodology.
    tds = torusinvariant_prime_divisors(ambient_space(m))
    cds = [cohomology_class(td) for td in tds]
    pt_class = cohomology_class(anticanonical_divisor_class(ambient_space(m)))
    for i in 1:length(ambient_space_flux_candidates_basis)
      condition = Vector{ZZRingElem}()
      for j in 1:length(list_of_divisor_pairs_to_be_considered)
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., list_of_divisor_pairs_to_be_considered[j]...]))
        push!(condition, get!(inter_dict, my_tuple) do
          class = ambient_space_flux_candidates_basis[i] * cds[list_of_divisor_pairs_to_be_considered[j][1]] * cds[list_of_divisor_pairs_to_be_considered[j][2]] * pt_class
          return ZZ(integrate(class))
        end)
      end
      push!(quant_constraint_matrix, condition)
    end

  end
  
  # (8) Convert the quant_constraint_matrix to a ZZ matrix. If necessary, multiply it by a suitable integer.
  denom = lcm(unique(vcat([denominator.(k) for k in quant_constraint_matrix]...)))
  if denom != 1
    quant_constraint_matrix = denom * quant_constraint_matrix
  end
  C = transpose(matrix(ZZ, quant_constraint_matrix))


  # (9) Work out the well-quantized fluxes as linear combinations of the parametrization of the vertical fluxes
  C2 = C * vertical_fluxes # This is a ZZ-matrix, since we parametrize vertical fluxes with integer coefficients!
  S, T, U = snf_with_transform(C2)
  r = rank(S)
  @req all(k -> !is_zero(S[k,k]), 1:r) "Inconsistency in Smith normal form detected. Please inform the authors."
  @req all(k -> is_zero(S[k,k]), r+1:min(nrows(S), ncols(S))) "Inconsistency in Smith normal form  detected. Please inform the authors."
  S_prime = zero_matrix(QQ, ncols(S), ncols(S))
  for k in 1:min(nrows(S), ncols(S))
    if k <= r
      S_prime[k,k] = denom//S[k,k]
    else
      S_prime[k,k] = 1
    end
  end
  solution_matrix = U * S_prime


  # (10) Finally, we need to re-express those in terms of the original bases.
  # (10) Rather, we have res now expressed in terms of the basis of vertical fluxes...
  sol_mat = vertical_fluxes * solution_matrix
  res = (sol_mat[:,1:r], sol_mat[:,r+1:ncols(solution_matrix)])


  # (11) Remember computed data
  fgs = family_of_g4_fluxes(m, res[1], res[2])
  set_attribute!(m, :matrix_integral_quant_transverse, res[1])
  set_attribute!(m, :matrix_rational_quant_transverse, res[2])
  set_attribute!(fgs, :is_well_quantized, true)
  set_attribute!(fgs, :passes_transversality_checks, true)
  set_attribute!(fgs, :breaks_non_abelian_gauge_group, true)
  set_attribute!(m, :inter_dict, inter_dict)
  set_attribute!(m, :s_inter_dict, s_inter_dict)


  # (12) Finally, return the result
  return fgs
end


@attr FamilyOfG4Fluxes function well_quantized_and_transversal_and_no_non_abelian_gauge_group_breaking_ambient_space_models_of_g4_fluxes(m::AbstractFTheoryModel; check::Bool = true)

  # (0) Has this result been computed before?
  if has_attribute(m, :matrix_integral_quant_transverse_nobreak) && has_attribute(m, :matrix_rational_quant_transverse_nobreak)
    fgs = family_of_g4_fluxes(m, matrix_integral_quant_transverse_nobreak(m, check = check), matrix_rational_quant_transverse_nobreak(m, check = check))
    set_attribute!(fgs, :is_well_quantized, true)
    set_attribute!(fgs, :passes_transversality_checks, true)
    set_attribute!(fgs, :breaks_non_abelian_gauge_group, false)
    return fgs
  end

  
  # (1) Entry checks
  @req base_space(m) isa NormalToricVariety "Computation of well-quantized, transversal and non-breaking G4-fluxes only supported for toric base and ambient spaces"
  @req dim(ambient_space(m)) == 5 "Computation of well-quantized, transversal and non-breaking G4-fluxes only supported for 5-dimensional toric ambient spaces"
  if check
    @req is_complete(ambient_space(m)) "Computation of well-quantized, transversal and non-breaking G4-fluxes only supported for complete toric ambient spaces"
    @req is_simplicial(ambient_space(m)) "Computation of well-quantized, transversal and non-breaking G4-fluxes only supported for simplicial toric ambient space"
  end


  # (2) Compute data, that is frequently used by the sophisticated intersection product below
  S = cox_ring(ambient_space(m))
  gS = gens(cox_ring(ambient_space(m)))
  linear_relations = matrix(QQ, rays(ambient_space(m)))
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
  # TODO: If available and necessary, convert inter_dict.
  # TODO: This is necessary, because serializing and loading turns NTuple{4, Int64} into Tuple (as of March 5, 2025).
  # TODO: Once serialization has caught up, this conversion will no longer be needed.
  if has_attribute(m, :inter_dict) && typeof(get_attribute(m, :inter_dict)) != Dict{NTuple{4, Int64}, ZZRingElem}
    original_dict = get_attribute(m, :inter_dict)
    new_dict = Dict{NTuple{4, Int64}, ZZRingElem}()
    for (key, value) in original_dict
      new_key = NTuple{4, Int64}(key)
      new_dict[new_key] = value
    end
    set_attribute!(model, :inter_dict, new_dict)
  end
  inter_dict = get_attribute!(m, :inter_dict) do
    Dict{NTuple{4, Int64}, ZZRingElem}()
  end::Dict{NTuple{4, Int64}, ZZRingElem}
  s_inter_dict = get_attribute!(m, :s_inter_dict) do
    Dict{String, ZZRingElem}()
  end::Dict{String, ZZRingElem}


  # (4) Obtain critical information - this may take significant time!
  ambient_space_flux_candidates_basis = _ambient_space_models_of_g4_fluxes(m, check = check)
  list_of_base_divisor_pairs_to_be_considered = Oscar._ambient_space_base_divisor_pairs_to_be_considered(m)
  ambient_space_flux_candidates_basis_indices = get_attribute(m, :ambient_space_models_of_g4_fluxes_indices)::Vector{Tuple{Int64, Int64}}
  list_of_divisor_pairs_to_be_considered = Oscar._ambient_space_divisor_pairs_to_be_considered(m)


  # (5) The following is fragile, but hopefully is a starting point
  exceptional_divisor_positions = findall(x -> occursin(r"^e\d+(_\d+)?$", x), string.(symbols(S)))


  # (6) Work out the relevant intersection numbers to tell if a flux is vertical
  vertical_and_no_gauge_group_breaking_constraint_matrix = Vector{Vector{QQFieldElem}}()
  # I have prepared some functionality below, regarding the case that this matrix should have rational entries.
  # However, I expect that this will not happen as long as the hypersurface in question is smooth.
  if arxiv_doi(m) == "10.48550/arXiv.1511.03209"

    # Use special intersection theory for special F-theory model. This technology could be extended beyond this one use-case in the future.
    for i in 1:length(ambient_space_flux_candidates_basis)

      condition = Vector{ZZRingElem}()

      # Compute against pairs of base divisors
      for j in 1:length(list_of_base_divisor_pairs_to_be_considered)
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., list_of_base_divisor_pairs_to_be_considered[j]...]))
        push!(condition, sophisticated_intersection_product(ambient_space(m), my_tuple, hypersurface_equation(m), inter_dict, s_inter_dict, data))
      end

      # Compute against zero section and base divisor
      pos_zero_section = findfirst(x -> x == "z", string.(gS))
      for j in 1:n_rays(base_space(m))
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., [j, pos_zero_section]...]))
        push!(condition, sophisticated_intersection_product(ambient_space(m), my_tuple, hypersurface_equation(m), inter_dict, s_inter_dict, data))
      end

      # Compute against exceptional divisors
      for j in 1:n_rays(base_space(m))
        for k in 1:length(exceptional_divisor_positions)
          my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., [j, exceptional_divisor_positions[k]]...]))
          push!(condition, sophisticated_intersection_product(ambient_space(m), my_tuple, hypersurface_equation(m), inter_dict, s_inter_dict, data))
        end
      end

      # Remember the computed condition
      push!(vertical_and_no_gauge_group_breaking_constraint_matrix, condition)

    end

  else
  
    # Cover all other case with generic, but potentially painfully slow methodology.
    tds = torusinvariant_prime_divisors(ambient_space(m))
    cds = [cohomology_class(td) for td in tds]
    pt_class = cohomology_class(anticanonical_divisor_class(ambient_space(m)))
    for i in 1:length(ambient_space_flux_candidates_basis)

      condition = Vector{ZZRingElem}()

      # Compute against pairs of base divisors
      for j in 1:length(list_of_base_divisor_pairs_to_be_considered)
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., list_of_base_divisor_pairs_to_be_considered[j]...]))
        push!(condition, get!(inter_dict, my_tuple) do
          class = ambient_space_flux_candidates_basis[i] * cds[list_of_base_divisor_pairs_to_be_considered[j][1]] * cds[list_of_base_divisor_pairs_to_be_considered[j][2]] * pt_class
          return ZZ(integrate(class))
        end)
      end

      # Compute against zero section and base divisor
      zsc = zero_section_class(m)
      pos_zero_section = findfirst(x -> x == string(polynomial(zsc)), string.([polynomial(x) for x in cds]))
      @req pos_zero_section !== nothing && pos_zero_section >= 1 "Could not establish position of the zero section"
      for j in 1:n_rays(base_space(m))
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., [j, pos_zero_section]...]))
        push!(condition, get!(inter_dict, my_tuple) do
          class = ambient_space_flux_candidates_basis[i] * cds[j] * zsc * pt_class
          return ZZ(integrate(class))
        end)

      end

      # Compute against exceptional divisors
      for j in 1:n_rays(base_space(m))
        for k in 1:length(exceptional_divisor_positions)
          my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., j, exceptional_divisor_positions[k]...]))
          push!(condition, get!(inter_dict, my_tuple) do
            class = ambient_space_flux_candidates_basis[i] * cds[j] * cds[exceptional_divisor_positions[k]] * pt_class
            return ZZ(integrate(class))
          end)
        end
      end

      # Remember the computed condition
      push!(vertical_and_no_gauge_group_breaking_constraint_matrix, condition)

    end

  end


  # (6) Compute the vertical fluxes as the kernel of the vertical_and_no_gauge_group_breaking_constraint_matrix.
  # (6) To later tell if those fluxes are properly quantized, we want to parametrize them with integer coefficient only.
  denom = lcm(unique(vcat([denominator.(k) for k in vertical_and_no_gauge_group_breaking_constraint_matrix]...)))
  if denom != 1
    vertical_and_no_gauge_group_breaking_constraint_matrix = denom * vertical_and_no_gauge_group_breaking_constraint_matrix
  end
  C_vertical_and_no_gauge_group_breaking = transpose(matrix(ZZ, vertical_and_no_gauge_group_breaking_constraint_matrix))
  vertical_and_no_gauge_group_breaking_fluxes = nullspace(C_vertical_and_no_gauge_group_breaking)[2]


  # (7) Work out the relevant intersection numbers to tell if a flux is well quantized
  quant_constraint_matrix = Vector{Vector{QQFieldElem}}()
  # I have prepared some functionality below, regarding the case that this matrix should have rational entries.
  # However, I expect that this will not happen as long as the hypersurface in question is smooth.
  if arxiv_doi(m) == "10.48550/arXiv.1511.03209"

    # Use special intersection theory for special F-theory model. This technology could be extended beyond this one use-case in the future.
    for i in 1:length(ambient_space_flux_candidates_basis)
      condition = Vector{ZZRingElem}()
      for j in 1:length(list_of_divisor_pairs_to_be_considered)
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., list_of_divisor_pairs_to_be_considered[j]...]))
        push!(condition, sophisticated_intersection_product(ambient_space(m), my_tuple, hypersurface_equation(m), inter_dict, s_inter_dict, data))
      end
      push!(quant_constraint_matrix, condition)
    end

  else
  
    # Cover all other case with generic, but potentially painfully slow methodology.
    tds = torusinvariant_prime_divisors(ambient_space(m))
    cds = [cohomology_class(td) for td in tds]
    pt_class = cohomology_class(anticanonical_divisor_class(ambient_space(m)))
    for i in 1:length(ambient_space_flux_candidates_basis)
      condition = Vector{ZZRingElem}()
      for j in 1:length(list_of_divisor_pairs_to_be_considered)
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., list_of_divisor_pairs_to_be_considered[j]...]))
        push!(condition, get!(inter_dict, my_tuple) do
          class = ambient_space_flux_candidates_basis[i] * cds[list_of_divisor_pairs_to_be_considered[j][1]] * cds[list_of_divisor_pairs_to_be_considered[j][2]] * pt_class
          return ZZ(integrate(class))
        end)
      end
      push!(quant_constraint_matrix, condition)
    end

  end
  
  # (8) Convert the quant_constraint_matrix to a ZZ matrix. If necessary, multiply it by a suitable integer.
  denom = lcm(unique(vcat([denominator.(k) for k in quant_constraint_matrix]...)))
  if denom != 1
    quant_constraint_matrix = denom * quant_constraint_matrix
  end
  C = transpose(matrix(ZZ, quant_constraint_matrix))


  # (9) Work out the well-quantized fluxes as linear combinations of the parametrization of the vertical and no non-abelian gauge group breaking fluxes
  C2 = C * vertical_and_no_gauge_group_breaking_fluxes # This is a ZZ-matrix, since we parametrize vertical fluxes with integer coefficients!
  S, T, U = snf_with_transform(C2)
  r = rank(S)
  @req all(k -> !is_zero(S[k,k]), 1:r) "Inconsistency in Smith normal form detected. Please inform the authors."
  @req all(k -> is_zero(S[k,k]), r+1:min(nrows(S), ncols(S))) "Inconsistency in Smith normal form  detected. Please inform the authors."
  S_prime = zero_matrix(QQ, ncols(S), ncols(S))
  for k in 1:min(nrows(S), ncols(S))
    if k <= r
      S_prime[k,k] = denom//S[k,k]
    else
      S_prime[k,k] = 1
    end
  end
  solution_matrix = U * S_prime


  # (10) Finally, we need to re-express those in terms of the original bases.
  # (10) Rather, we have res now expressed in terms of the basis of vertical and no non-abelian gauge group breaking fluxes...
  sol_mat = vertical_and_no_gauge_group_breaking_fluxes * solution_matrix
  res = (sol_mat[:,1:r], sol_mat[:,r+1:ncols(solution_matrix)])


  # (11) Remember computed data
  fgs = family_of_g4_fluxes(m, res[1], res[2])
  set_attribute!(m, :matrix_integral_quant_transverse_nobreak, res[1])
  set_attribute!(m, :matrix_rational_quant_transverse_nobreak, res[2])
  set_attribute!(fgs, :is_well_quantized, true)
  set_attribute!(fgs, :passes_transversality_checks, true)
  set_attribute!(fgs, :breaks_non_abelian_gauge_group, false)
  set_attribute!(m, :inter_dict, inter_dict)
  set_attribute!(m, :s_inter_dict, s_inter_dict)


  # (12) Finally, return the result
  return fgs
end
