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

  # Each basis element is given by the vanishing of two homogeneous variables. We extract those indices...
  filtered_h22_basis_indices_init = get_attribute(ambient_space(m), :basis_of_h22_indices)

  # It may happen that filtered_h22_basis_indices_init is encoded as Vector{Any}. But it is a Vector{Tuple{Int64, Int64}}
  # Of course, this should be fixed more properly, but for now, the following works...
  filtered_h22_basis_indices = [k for k in filtered_h22_basis_indices_init]::Vector{Tuple{Int64, Int64}}

  # Prepare data of the toric ambient space
  gS = gens(cox_ring(ambient_space(m)))
  mnf = Oscar._minimal_nonfaces(ambient_space(m))
  sr_ideal_pos = [Vector{Int}(Polymake.row(mnf, i)) for i in 1:Polymake.nrows(mnf)]

  # Filter out basis elements
  for a in length(filtered_h22_basis):-1:1
    
    # Simplify the hypersurface polynomial by setting relevant variables to zero
    vanishing_vars_pos = [filtered_h22_basis_indices[a]...]
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
      deleteat!(filtered_h22_basis_indices, a)
    end

  end

  set_attribute!(m, :ambient_space_models_of_g4_fluxes, filtered_h22_basis)
  set_attribute!(m, :ambient_space_models_of_g4_fluxes_indices, filtered_h22_basis_indices)
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

Crucially, this method assumes that $c_2( \widehat{Y}_4)$ is even. Currently, no
check is conducted and no error raised. Use with care!

Recall that this is relevant in so much as the quantization condition asks to
verify if the twist of the given $G_4$-flux by $1/2 \cdot c_2( \widehat{Y}_4)$ is
even. Recall also that it is known that for many F-theory models, $c_2( \widehat{Y}_4)$
is an even class. For instance, this applies to all F-theory compactifications
on an elliptically fibered smooth Calabi-Yau 4-fold with a globally defined
Weierstrass model [CS12](@cite). For instance, this means that all of the
F-theory QSMs [CHLLT19](@cite) have an even $c_2( \widehat{Y}_4)$.

It can be computationally very demanding to check if a toric variety
$X$ is complete (and simplicial). The optional argument `check` can be set
to `false` to skip these tests.

# Examples
```jldoctest
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

Here is a more interesting example.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> g4_base = ambient_space_models_of_g4_fluxes(qsm_model, check = false);

julia> length(g4_base)
37

julia> res = well_quantized_ambient_space_models_of_g4_fluxes(qsm_model, check = false);

julia> size(res[1])
(37, 37)

julia> size(res[2])
(37, 0)

julia> M = res[1];

julia> g4_class = sum(M[i,j]*g4_base[i] for i in 1:length(g4_base) for j in 1:size(M,2));

julia> g4 = g4_flux(qsm_model, g4_class, check = false)
G4-flux candidate lacking elementary quantization checks

julia> passes_elementary_quantization_checks(g4)
true
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
  inter_dict = Dict{NTuple{4, Int64}, ZZRingElem}()
  s_inter_dict = Dict{String, ZZRingElem}()
  if has_attribute(m, :inter_dict)
    inter_dict = get_attribute(m, :inter_dict)
  end
  if has_attribute(m, :s_inter_dict)
    s_inter_dict = get_attribute(m, :s_inter_dict)
  end


  # (4) Obtain critical information - this may take significant time!
  ambient_space_flux_candidates_basis = ambient_space_models_of_g4_fluxes(m, check = check)
  ambient_space_flux_candidates_basis_indices = get_attribute(m, :ambient_space_models_of_g4_fluxes_indices)
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
    cds = [cohomology_class(tds[i]) for i in 1:length(tds)]
    pt_class = cohomology_class(anticanonical_divisor_class(ambient_space(m)))
    for i in 1:length(ambient_space_flux_candidates_basis)
      condition = Vector{ZZRingElem}()
      for j in 1:length(list_of_divisor_pairs_to_be_considered)
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., list_of_divisor_pairs_to_be_considered[j]...]))
        if !haskey(inter_dict, my_tuple)
          class = ambient_space_flux_candidates_basis[i] * cds[list_of_divisor_pairs_to_be_considered[j][1]] * cds[list_of_divisor_pairs_to_be_considered[j][2]] * pt_class
          inter_dict[my_tuple] = ZZ(integrate(class))
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


@doc raw"""
    well_quantized_and_vertical_ambient_space_models_of_g4_fluxes(m::AbstractFTheoryModel; check::Bool = true)::Tuple{QQMatrix, QQMatrix}

Given an F-theory model $m$ defined as hypersurface in a simplicial and
complete toric base, this method computes a basis of all well-quantized
and vertical  ambient space $G_4$-fluxes. The result of this operation is a
tuple of two matrices. The columns of the first matrix specify those (rational)
combinations of ambient space $G_4$-fluxes, of which one may only take
$\mathbb{Z}$-linear combinations without violating flux quantization. The columns
of the second matrix specify those (rational) combinations of ambient space
$G_4$-fluxes, for which any rational linear combination satisfies the flux
quantization condition.

Crucially, this method assumes that $c_2( \widehat{Y}_4)$ is even. Currently, no
check is conducted and no error raised. Use with care!

Recall that this is relevant in so much as the quantization condition asks to
verify if the twist of the given $G_4$-flux by $1/2 \cdot c_2( \widehat{Y}_4)$ is
even. Recall also that it is known that for many F-theory models, $c_2( \widehat{Y}_4)$
is an even class. For instance, this applies to all F-theory compactifications
on an elliptically fibered smooth Calabi-Yau 4-fold with a globally defined
Weierstrass model [CS12](@cite). For instance, this means that all of the
F-theory QSMs [CHLLT19](@cite) have an even $c_2( \widehat{Y}_4)$.

It can be computationally very demanding to check if a toric variety
$X$ is complete (and simplicial). The optional argument `check` can be set
to `false` to skip these tests.

# Examples
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> Kbar = anticanonical_divisor_class(B3)
Divisor class on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w"=>Kbar))
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> res = well_quantized_and_vertical_ambient_space_models_of_g4_fluxes(t, check = false);

julia> res[1]
2 by 0 empty matrix

julia> res[2]
2 by 0 empty matrix
```

Here is a more interesting example, in which we verify with our software tool for one particular F-theory QSM, that the
choice of $G_4$-flux presented in [CHLLT19](@cite), is indeed vertical and satisfies the necessary conditions
for being well-quantized.

```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 2021))
Hypersurface model over a concrete base

julia> res = well_quantized_and_vertical_ambient_space_models_of_g4_fluxes(qsm_model, check = false);

julia> size(res[1])
(37, 25)

julia> size(res[2])
(37, 0)

julia> M=res[1];

julia> g4_base = ambient_space_models_of_g4_fluxes(qsm_model, check = false);

julia> g4_classes = [sum(M[i,j]*g4_base[i] for i in 1:length(g4_base)) for j in 1:size(M,2)];

julia> length(g4_classes) == 25
true

julia> g4_classes[end]
Cohomology class on a normal toric variety given by 293//300*x4*e2 + 143//150*x4*u - 283//25*x4*e4 + 143//150*x4*e1 + 1643//300*x4*w - 599//150*x5*x8 - 7//150*x5*e2 - 7//75*x5*u - 7//50*x5*e4 - 7//75*x5*e1 - 89//300*x5*w + 896//75*x6*x7 + 1//20*x6*e2 + 1//10*x6*u + 2//5*x6*e4 + 1//10*x6*e1 - 1//5*x6*w - 599//75*x7*x8 - 7//75*x7*e2 - 14//75*x7*u - 7//25*x7*e4 - 14//75*x7*e1 - 89//150*x7*w + 208//75*x8^2 + 298//75*x8*x9 + 1//150*x8*e2 - 73//150*x8*u - 12//25*x8*e4 - 73//150*x8*e1 - 37//75*x8*w + 82//75*x9^2 - 7//150*x9*e2 - 7//75*x9*u + 9//25*x9*e4 - 7//75*x9*e1 - 41//75*x9*w + 11//30*e1*w

julia> g4_list = [g4_flux(qsm_model, cl, check = false) for cl in g4_classes];

julia> all(k -> passes_elementary_quantization_checks(k), g4_list)
true

julia> all(k -> passes_verticality_checks(k), g4_list)
true

julia> c = [60, 51, 90, 0, 24, 51, -24, 45, 30, 0, -48, 90, -57, 60, 30, 15, 120, 0, -60, 0, -720, -420, -270, -60, -2190];

julia> qsm_g4_candidate = g4_flux(qsm_model, sum(c[i]*g4_classes[i] for i in 1:length(g4_classes)), check = false)
G4-flux candidate lacking elementary quantization checks

julia> passes_elementary_quantization_checks(qsm_g4_candidate)
true

julia> passes_verticality_checks(qsm_g4_candidate)
true
```
"""
function well_quantized_and_vertical_ambient_space_models_of_g4_fluxes(m::AbstractFTheoryModel; check::Bool = true)::Tuple{QQMatrix, QQMatrix}

  # (1) Entry checks
  @req base_space(m) isa NormalToricVariety "Computation of well-quantized G4-fluxes only supported for toric base and ambient spaces"
  @req dim(ambient_space(m)) == 5 "Computation of well-quantized G4-fluxes only supported for 5-dimensional toric ambient spaces"
  if check
    @req is_complete(ambient_space(m)) "Computation of well-quantized G4-fluxes only supported for complete toric ambient spaces"
    @req is_simplicial(ambient_space(m)) "Computation of well-quantized G4-fluxes only supported for simplicial toric ambient space"
  end
  if has_attribute(m, :well_quantized_and_vertical_ambient_space_models_of_g4_fluxes)
    return get_attribute(m, :well_quantized_and_vertical_ambient_space_models_of_g4_fluxes)
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
  inter_dict = Dict{NTuple{4, Int64}, ZZRingElem}()
  s_inter_dict = Dict{String, ZZRingElem}()
  if has_attribute(m, :inter_dict)
    inter_dict = get_attribute(m, :inter_dict)
  end
  if has_attribute(m, :s_inter_dict)
    s_inter_dict = get_attribute(m, :s_inter_dict)
  end


  # (4) Obtain critical information - this may take significant time!
  ambient_space_flux_candidates_basis = ambient_space_models_of_g4_fluxes(m, check = check)
  list_of_base_divisor_pairs_to_be_considered = Oscar._ambient_space_base_divisor_pairs_to_be_considered(m)
  ambient_space_flux_candidates_basis_indices = get_attribute(m, :ambient_space_models_of_g4_fluxes_indices)
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
    cds = [cohomology_class(tds[i]) for i in 1:length(tds)]
    pt_class = cohomology_class(anticanonical_divisor_class(ambient_space(m)))
    for i in 1:length(ambient_space_flux_candidates_basis)

      condition = Vector{ZZRingElem}()

      # Compute against pairs of base divisors
      for j in 1:length(list_of_base_divisor_pairs_to_be_considered)
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., list_of_base_divisor_pairs_to_be_considered[j]...]))
        if !haskey(inter_dict, my_tuple)
          class = ambient_space_flux_candidates_basis[i] * cds[list_of_base_divisor_pairs_to_be_considered[j][1]] * cds[list_of_base_divisor_pairs_to_be_considered[j][2]] * pt_class
          inter_dict[my_tuple] = ZZ(integrate(class))
        end
        push!(condition, inter_dict[my_tuple])
      end

      # Compute against zero section and base divisor
      zsc = zero_section_class(m)
      pos_zero_section = findfirst(x -> x == string(polynomial(zsc)), string.([polynomial(x) for x in cds]))
      @req pos_zero_section !== nothing && pos_zero_section >= 1 "Could not establish position of the zero section"
      for j in 1:n_rays(base_space(m))
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., [j, pos_zero_section]...]))
        if !haskey(inter_dict, my_tuple)  
          class = ambient_space_flux_candidates_basis[i] * cds[j] * zsc * pt_class
          inter_dict[my_tuple] = ZZ(integrate(class))
        end
        push!(condition, inter_dict[my_tuple])

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
    cds = [cohomology_class(tds[i]) for i in 1:length(tds)]
    pt_class = cohomology_class(anticanonical_divisor_class(ambient_space(m)))
    for i in 1:length(ambient_space_flux_candidates_basis)
      condition = Vector{ZZRingElem}()
      for j in 1:length(list_of_divisor_pairs_to_be_considered)
        my_tuple = Tuple(sort([ambient_space_flux_candidates_basis_indices[i]..., list_of_divisor_pairs_to_be_considered[j]...]))
        if !haskey(inter_dict, my_tuple)
          class = ambient_space_flux_candidates_basis[i] * cds[list_of_divisor_pairs_to_be_considered[j][1]] * cds[list_of_divisor_pairs_to_be_considered[j][2]] * pt_class
          inter_dict[my_tuple] = ZZ(integrate(class))
        end
        push!(condition, inter_dict[my_tuple])
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
  set_attribute!(m, :well_quantized_and_vertical_ambient_space_models_of_g4_fluxes, res)
  set_attribute!(m, :inter_dict, inter_dict)
  set_attribute!(m, :s_inter_dict, s_inter_dict)


  # (12) Finally, return the result
  return res

end
