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

julia> qsm_g4_flux == qsm_flux(qsm_model)
true
```
"""
function special_flux_family(m::AbstractFTheoryModel; not_breaking::Bool = false, check::Bool = true, algorithm::String = "default")

  # (1) Consistency checks
  if arxiv_doi(m) == "10.48550/arXiv.1511.03209" && algorithm == "default"
    error("The default algorithm for intersection computations will likely not terminate in a reasonable time for this model and is therefore not supported")
  end
  @req base_space(m) isa NormalToricVariety "Computation of well-quantized and transversal G4-fluxes only supported for toric base and ambient spaces"
  @req dim(ambient_space(m)) == 5 "Computation of well-quantized and transversal G4-fluxes only supported for 5-dimensional toric ambient spaces"
  if check
    @req is_complete(ambient_space(m)) "Computation of well-quantized and transversal G4-fluxes only supported for complete toric ambient spaces"
    @req is_simplicial(ambient_space(m)) "Computation of well-quantized and transversal G4-fluxes only supported for simplicial toric ambient space"
  end

  # (2) Is result known?
  if !not_breaking
    if has_attribute(m, :matrix_integral_quant_transverse) && has_attribute(m, :matrix_rational_quant_transverse) && has_attribute(m, :offset_quant_transverse)
      fgs_m_int = matrix_integral_quant_transverse(m, check = check)
      fgs_m_rat = matrix_rational_quant_transverse(m, check = check)
      fgs_offset = offset_quant_transverse(m, check = check)
      fgs = family_of_g4_fluxes(m, fgs_m_int, fgs_m_rat, fgs_offset)
      set_attribute!(fgs, :is_well_quantized, true)
      set_attribute!(fgs, :passes_transversality_checks, true)
      set_attribute!(fgs, :breaks_non_abelian_gauge_group, true)
      return fgs
    end
  else
    if has_attribute(m, :matrix_integral_quant_transverse_nobreak) && has_attribute(m, :matrix_rational_quant_transverse_nobreak) && has_attribute(m, :offset_quant_transverse_nobreak)
      fgs_m_int = matrix_integral_quant_transverse_nobreak(m, check = check)
      fgs_m_rat = matrix_rational_quant_transverse_nobreak(m, check = check)
      fgs_offset = offset_quant_transverse_nobreak(m, check = check)
      fgs = family_of_g4_fluxes(m, fgs_m_int, fgs_m_rat, fgs_offset)
      set_attribute!(fgs, :is_well_quantized, true)
      set_attribute!(fgs, :passes_transversality_checks, true)
      set_attribute!(fgs, :breaks_non_abelian_gauge_group, false)
      return fgs
    end  
  end

  # (3) Result not known, compute it!
  final_shift = Vector{QQFieldElem}()
  res = Vector{ZZMatrix}()
  if algorithm == "special"
    final_shift, res = special_flux_family_with_special_algorithm(m::AbstractFTheoryModel; not_breaking = not_breaking, check = check)
  else
    final_shift, res = special_flux_family_with_default_algorithm(m::AbstractFTheoryModel; not_breaking = not_breaking, check = check)
  end

  # (4) Set attributes accordingly, and return result
  fgs = family_of_g4_fluxes(m, res[1], res[2], final_shift)
  set_attribute!(fgs, :is_well_quantized, true)
  set_attribute!(fgs, :passes_transversality_checks, true)
  if !not_breaking
    set_attribute!(m, :matrix_integral_quant_transverse, res[1])
    set_attribute!(m, :matrix_rational_quant_transverse, res[2])
    set_attribute!(m, :offset_quant_transverse, final_shift)
    set_attribute!(fgs, :breaks_non_abelian_gauge_group, true)
  else
    set_attribute!(m, :matrix_integral_quant_transverse_nobreak, res[1])
    set_attribute!(m, :matrix_rational_quant_transverse_nobreak, res[2])
    set_attribute!(m, :offset_quant_transverse_nobreak, final_shift)
    set_attribute!(fgs, :breaks_non_abelian_gauge_group, false)
  end
  return fgs

end


function special_flux_family_with_default_algorithm(m::AbstractFTheoryModel; not_breaking::Bool = false, check::Bool = true)  
  
  # (1) Are intersection numbers known?
  # These instructions appear twice, once in the default and once in the special algorithnm. Code duplication? Improve it!
  inter_dict = get_attribute!(m, :inter_dict) do
    Dict{NTuple{4, Int64}, ZZRingElem}()
  end::Dict{NTuple{4, Int64}, ZZRingElem}
  s_inter_dict = get_attribute!(m, :s_inter_dict) do
    Dict{String, ZZRingElem}()
  end::Dict{String, ZZRingElem}


  # (2) Obtain critical information - this may take significant time!
  ambient_space_flux_candidates_basis = basis_of_h22_hypersurface(m, check = check)
  list_of_base_divisor_pairs_to_be_considered = Oscar._ambient_space_base_divisor_pairs_to_be_considered(m)
  ambient_space_flux_candidates_basis_indices = basis_of_h22_hypersurface_indices(m, check = check)
  list_of_divisor_pairs_to_be_considered = Oscar._ambient_space_divisor_pairs_to_be_considered(m)
  S = cox_ring(ambient_space(m))
  exceptional_divisor_positions = findall(x -> occursin(r"^e\d+(_\d+)?$", x), string.(symbols(S))) # TODO: This line is a bit fragile. Fix it!
  tds = torusinvariant_prime_divisors(ambient_space(m))
  cds = [cohomology_class(td) for td in tds]
  pt_class = cohomology_class(anticanonical_divisor_class(ambient_space(m)))


  # (3) Work out the relevant intersection numbers to tell if a flux passes the transversality constraints & (if desired) does not break the non-abelian gauge group.
  transversality_constraint_matrix = Vector{Vector{ZZRingElem}}()
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
    pos_zero_section = zero_section_index(m)
    @req pos_zero_section !== nothing && pos_zero_section >= 1 "Could not establish position of the zero section"
    for j in 1:n_rays(base_space(m))
      my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., [j, pos_zero_section]...]))
      push!(condition, get!(inter_dict, my_tuple) do
        class = ambient_space_flux_candidates_basis[i] * cds[j] * zsc * pt_class
        return ZZ(integrate(class))
      end)
    end

    # Compute against exceptional divisors
    if not_breaking
      for j in 1:n_rays(base_space(m))
        for k in 1:length(exceptional_divisor_positions)
          my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., j, exceptional_divisor_positions[k]...]))
          push!(condition, get!(inter_dict, my_tuple) do
            class = ambient_space_flux_candidates_basis[i] * cds[j] * cds[exceptional_divisor_positions[k]] * pt_class
            return ZZ(integrate(class))
          end)
        end
      end
    end      

    # Remember the computed result
    push!(transversality_constraint_matrix, condition)

  end
  C_transverse = transpose(matrix(ZZ, transversality_constraint_matrix))
  transverse_fluxes = nullspace(C_transverse)[2]


  # (4) Work out the relevant intersection numbers to tell if a flux is well quantized
  quant_constraint_matrix = Vector{Vector{ZZRingElem}}()
  offset_vector = Vector{QQFieldElem}()
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
  C = transpose(matrix(ZZ, quant_constraint_matrix))
  # TODO: We currently do not remember any intersection numbers computed from this operation. Improve.
  offset = 1//2 * chern_class(m, 2, check = check)
  for j in 1:length(list_of_divisor_pairs_to_be_considered)
    class = offset * cds[list_of_divisor_pairs_to_be_considered[j][1]] * cds[list_of_divisor_pairs_to_be_considered[j][2]] * pt_class
    push!(offset_vector, QQ(integrate(class)))
  end
  

  # (5) Work out the well-quantized fluxes as linear combinations of the parametrization of the fluxes which pass the transversality constraints.
  C2 = C * transverse_fluxes # Intersection numbers in terms off the basis of transverse fluxes.
  S, T, U = snf_with_transform(C2)
  r = rank(S)
  @req all(k -> !is_zero(S[k,k]), 1:r) "Inconsistency in Smith normal form computation detected. Please inform the authors."
  @req all(k -> is_zero(S[k,k]), r+1:min(nrows(S), ncols(S))) "Inconsistency in Smith normal form computation detected. Please inform the authors."
  S_prime = zero_matrix(QQ, ncols(S), ncols(S))
  transformed_offset_vector = T * offset_vector
  shift_vector = [zero(QQ) for k in 1:nrows(S_prime)]
  for k in 1:min(nrows(S), ncols(S))
    if k <= r
      S_prime[k,k] = 1//S[k,k]
      if !isinteger(transformed_offset_vector[k])
        shift_vector[k] = - transformed_offset_vector[k]//S[k,k]
      end
    else
      @req isinteger(transformed_offset_vector[k]) "Inconsistency in Smith normal form computation detected. Please inform the authors."
      S_prime[k,k] = 1
    end
  end
  solution_matrix = U * S_prime
  solution_shift = U * shift_vector
  sol_mat = transverse_fluxes * solution_matrix
  final_shift = transverse_fluxes * solution_shift
  res = [sol_mat[:,1:r], sol_mat[:,r+1:ncols(solution_matrix)]]
  return [final_shift, res]

end


function special_flux_family_with_special_algorithm(m::AbstractFTheoryModel; not_breaking::Bool = false, check::Bool = true)
  
  # (1) Compute data, that is frequently used by the sophisticated intersection product below
  S = cox_ring(ambient_space(m))
  gS = gens(cox_ring(ambient_space(m)))
  linear_relations = matrix(ZZ, rays(ambient_space(m)))
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


  # (2) Are intersection numbers known?
  inter_dict = get_attribute!(m, :inter_dict) do
    Dict{NTuple{4, Int64}, ZZRingElem}()
  end::Dict{NTuple{4, Int64}, ZZRingElem}
  s_inter_dict = get_attribute!(m, :s_inter_dict) do
    Dict{String, ZZRingElem}()
  end::Dict{String, ZZRingElem}


  # (4) Obtain critical information - this may take significant time!
  ambient_space_flux_candidates_basis = basis_of_h22_hypersurface(m, check = check)
  list_of_base_divisor_pairs_to_be_considered = Oscar._ambient_space_base_divisor_pairs_to_be_considered(m)
  ambient_space_flux_candidates_basis_indices = basis_of_h22_hypersurface_indices(m, check = check)
  list_of_divisor_pairs_to_be_considered = Oscar._ambient_space_divisor_pairs_to_be_considered(m)
   # TODO: This line is a bit fragile. Fix it!
  exceptional_divisor_positions = findall(x -> occursin(r"^e\d+(_\d+)?$", x), string.(symbols(S)))


  # (5) Work out the relevant intersection numbers to tell if a flux passes the transversality constraints & (if desired) if the flux is not breaking the gauge group.
  transversality_constraint_matrix = Vector{Vector{ZZRingElem}}()
  for i in 1:length(ambient_space_flux_candidates_basis)

    condition = Vector{ZZRingElem}()

    # Compute against pairs of base divisors
    for j in 1:length(list_of_base_divisor_pairs_to_be_considered)
      my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., list_of_base_divisor_pairs_to_be_considered[j]...]))
      push!(condition, sophisticated_intersection_product(ambient_space(m), my_tuple, hypersurface_equation(m), inter_dict, s_inter_dict, data))
    end

    # Compute against zero section and base divisor
    pos_zero_section = zero_section_index(m)
    for j in 1:n_rays(base_space(m))
      my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., [j, pos_zero_section]...]))
      push!(condition, sophisticated_intersection_product(ambient_space(m), my_tuple, hypersurface_equation(m), inter_dict, s_inter_dict, data))
    end

    # Compute against exceptional divisors if desired
    if not_breaking
      for j in 1:n_rays(base_space(m))
        for k in 1:length(exceptional_divisor_positions)
          my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., [j, exceptional_divisor_positions[k]]...]))
          push!(condition, sophisticated_intersection_product(ambient_space(m), my_tuple, hypersurface_equation(m), inter_dict, s_inter_dict, data))
        end
      end
    end

    # Remember the computed intersection numbers
    push!(transversality_constraint_matrix, condition)

  end
  C_transverse = transpose(matrix(ZZ, transversality_constraint_matrix))
  transverse_fluxes = nullspace(C_transverse)[2]


  # (6) Work out the relevant intersection numbers to tell if a flux is well quantized
  quant_constraint_matrix = Vector{Vector{ZZRingElem}}()
  offset_vector = Vector{QQFieldElem}()
  for i in 1:length(ambient_space_flux_candidates_basis)
    condition = Vector{ZZRingElem}()
    for j in 1:length(list_of_divisor_pairs_to_be_considered)
      my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., list_of_divisor_pairs_to_be_considered[j]...]))
      push!(condition, sophisticated_intersection_product(ambient_space(m), my_tuple, hypersurface_equation(m), inter_dict, s_inter_dict, data))
    end
    push!(quant_constraint_matrix, condition)
  end
  C = transpose(matrix(ZZ, quant_constraint_matrix))
  c2 = lift(polynomial(chern_class(m, 2, check = check)))
  coeffs = collect(coefficients(c2))
  M = collect(exponents(lift(c2)))
  non_zero_exponents = Vector{Tuple{Int64, Int64}}()
  for my_row in M
    i1 = findfirst(x -> x != 0, my_row)
    my_row[i1] -= 1
    i2 = findfirst(x -> x != 0, my_row)
    push!(non_zero_exponents, (i1, i2))
  end
  for j in 1:length(list_of_divisor_pairs_to_be_considered)
    inter_numb = QQ(0)
    for k in 1:length(non_zero_exponents)
      my_tuple = Tuple(sort([non_zero_exponents[k]..., list_of_divisor_pairs_to_be_considered[j]...]))
      inter_numb += coeffs[k] * sophisticated_intersection_product(ambient_space(m), my_tuple, hypersurface_equation(m), inter_dict, s_inter_dict, data)
    end
    push!(offset_vector, inter_numb)
  end

  
  # (7) Work out the well-quantized fluxes as linear combinations of the parametrization of the fluxes which pass the transversality constraints.
  C2 = C * transverse_fluxes # Intersection numbers in terms off the basis of transverse fluxes.
  S, T, U = snf_with_transform(C2)
  r = rank(S)
  @req all(k -> !is_zero(S[k,k]), 1:r) "Inconsistency in Smith normal form computation detected. Please inform the authors."
  @req all(k -> is_zero(S[k,k]), r+1:min(nrows(S), ncols(S))) "Inconsistency in Smith normal form computation detected. Please inform the authors."
  S_prime = zero_matrix(QQ, ncols(S), ncols(S))
  transformed_offset_vector = T * offset_vector
  shift_vector = [zero(QQ) for k in 1:nrows(S_prime)]
  for k in 1:min(nrows(S), ncols(S))
    if k <= r
      S_prime[k,k] = 1//S[k,k]
      if !isinteger(transformed_offset_vector[k])
        shift_vector[k] = - transformed_offset_vector[k]//S[k,k]
      end
    else
      @req isinteger(transformed_offset_vector[k]) "Inconsistency in Smith normal form computation detected. Please inform the authors."
      S_prime[k,k] = 1
    end
  end
  solution_matrix = U * S_prime
  solution_shift = U * shift_vector
  sol_mat = transverse_fluxes * solution_matrix
  final_shift = transverse_fluxes * solution_shift
  res = [sol_mat[:,1:r], sol_mat[:,r+1:ncols(solution_matrix)]]
  return [final_shift, res]

end
