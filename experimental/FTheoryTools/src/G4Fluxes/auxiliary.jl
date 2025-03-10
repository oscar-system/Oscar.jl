@doc raw"""
    basis_of_h22(v::NormalToricVariety; check::Bool = true)::Vector{CohomologyClass}

By virtue of Theorem 12.4.1 in [CLS11](@cite), one can compute a monomial
basis of $H^4(X, \mathbb{Q})$ for a simplicial, complete toric variety $X$
by truncating its cohomology ring to degree $2$. Inspired by this, this
method identifies a basis of $H^{(2,2)}(X, \mathbb{Q})$ by multiplying
pairs of cohomology classes associated with toric coordinates.

By definition, $H^{(2,2)}(X, \mathbb{Q})$ is a subset of $H^{4}(X, \mathbb{Q})$.
However, by Theorem 9.3.2 in [CLS11](@cite), for complete and simplicial
toric varieties and $p \neq q$ it holds $H^{(p,q)}(X, \mathbb{Q}) = 0$. It follows
that for such varieties $H^{(2,2)}(X, \mathbb{Q}) = H^4(X, \mathbb{Q})$ and the
vector space dimension of those spaces agrees with the Betti number $b_4(X)$.

Note that it can be computationally very demanding to check if a toric variety
$X$ is complete (and simplicial). The optional argument `check` can be set
to `false` to skip these tests.

# Examples
```jldoctest
julia> Y1 = hirzebruch_surface(NormalToricVariety, 2)
Normal toric variety

julia> Y2 = hirzebruch_surface(NormalToricVariety, 2)
Normal toric variety

julia> Y = Y1 * Y2
Normal toric variety

julia> h22_basis = basis_of_h22(Y, check = false)
6-element Vector{CohomologyClass}:
 Cohomology class on a normal toric variety given by xx2*yx2
 Cohomology class on a normal toric variety given by xt2*yt2
 Cohomology class on a normal toric variety given by xx2*yt2
 Cohomology class on a normal toric variety given by xt2*yx2
 Cohomology class on a normal toric variety given by yx2^2
 Cohomology class on a normal toric variety given by xx2^2

julia> betti_number(Y, 4) == length(h22_basis)
true
```
"""
@attr Vector{CohomologyClass} function basis_of_h22(v::NormalToricVariety; check::Bool = true)

  # (0) Some initial checks
  if check
    @req is_complete(v) "Computation of basis of H22 is currently only supported for complete toric varieties"
    @req is_simplicial(v) "Computation of basis of H22 is currently only supported for simplicial toric varieties"
  end
  if dim(v) < 4
    return Vector{CohomologyClass}()
  end

  # (1) Prepare some data of the variety
  mnf = Oscar._minimal_nonfaces(v)
  ignored_sets = Set([Tuple(sort(Vector{Int}(Polymake.row(mnf, i)))) for i in 1:Polymake.nrows(mnf)])

  # (2) Prepare a dict that converts the "naive" generating set into our chosen basis
  converter_dict = Dict{Tuple{Int, Int}, Any}()
  for k in 1:n_rays(v)
    for l in k:n_rays(v)
      my_tuple = (k, l)
      if (my_tuple in ignored_sets)
        converter_dict[my_tuple] = 0
      else
        converter_dict[my_tuple] = nothing
      end
    end
  end

  # (3) Prepare the linear relations
  N_lin_rel, my_mat = rref(transpose(matrix(QQ, rays(v))))
  @req N_lin_rel == nrows(my_mat) "Cannot remove as many variables as there are linear relations - weird!"
  bad_positions = [findfirst(!iszero, row) for row in eachrow(my_mat)]
  lin_rels = Dict{Int, Vector{QQFieldElem}}()
  for k in 1:nrows(my_mat)
    my_relation = (-1) * my_mat[k, :]
    my_relation[bad_positions[k]] = 0
    @req all(k -> k == 0, my_relation[bad_positions]) "Inconsistency!"
    lin_rels[bad_positions[k]] = my_relation
  end

  # (4) Apply linear relations to remaining entries in converter_dict
  # (4) The code within the following for loop might need optimizing.
  crg = gens(base_ring(cohomology_ring(v)))
  for (key, value) in converter_dict
    if value === nothing
      image = Vector{Any}()

      # Only the first entry needs replacing by the linear relations
      if (key[1] in keys(lin_rels)) && (key[2] in keys(lin_rels)) == false
        my_relation = lin_rels[key[1]]
        posi = findall(k -> k != 0, my_relation)
        coeffs = my_relation[posi]
        for k in 1:length(posi)
          push!(image, [coeffs[k], (posi[k], key[2])])
        end
      end

      # Only the second entry needs replacing by the linear relations
      if (key[2] in keys(lin_rels)) && (key[1] in keys(lin_rels)) == false
        my_relation = lin_rels[key[2]]
        posi = findall(k -> k != 0, my_relation)
        coeffs = my_relation[posi]
        for k in 1:length(posi)
          push!(image, [coeffs[k], (key[1], posi[k])])
        end
      end

      # Both entry needs replacing by the linear relations
      if (key[1] in keys(lin_rels)) && key[2] in keys(lin_rels)
        my_relation = lin_rels[key[1]]
        posi = findall(k -> k != 0, my_relation)
        coeffs = my_relation[posi]
        for k in 1:length(posi)
          push!(image, [coeffs[k], (posi[k], key[2])])
        end
        my_relation = lin_rels[key[2]]
        posi = findall(k -> k != 0, my_relation)
        coeffs = my_relation[posi]
        old_image = copy(image)
        image = Vector{Any}()
        for i in 1:length(old_image)
          for k in 1:length(posi)
            old_coeff = old_image[i][1]
            push!(image, [old_coeff * coeffs[k], (old_image[i][2][1], posi[k])])
          end
        end
      end

      # If there was a replacement, update the key in the dict
      if length(image) > 0
        #image_as_cohomology_class = CohomologyClass(v, cohomology_ring(v)(sum(i[1] * crg[i[2][1]] * crg[i[2][2]] for i in image)))
        #converter_dict[key] = image_as_cohomology_class
        converter_dict[key] = image
      end
    end
  end

  # (5) Prepare a list of those variables that we keep, a.k.a. a basis of H^(1,1)
  good_positions = setdiff(1:n_rays(v), bad_positions)
  n_good_positions = length(good_positions)

  # (6) Make a list of all quadratic elements in the cohomology ring, which are not generators of the SR-ideal.
  N_filtered_quadratic_elements = 0
  dict_of_filtered_quadratic_elements = Dict{Tuple{Int64, Int64}, Int64}()
  for k in 1:n_good_positions
    for l in k:n_good_positions
      my_tuple = (min(good_positions[k], good_positions[l]), max(good_positions[k], good_positions[l]))
      if !(my_tuple in ignored_sets)
        N_filtered_quadratic_elements += 1
        dict_of_filtered_quadratic_elements[my_tuple] = N_filtered_quadratic_elements
      end
    end
  end

  # (7) Consistency check
  l1 = sort(collect(keys(converter_dict))[findall(k -> converter_dict[k] === nothing, collect(keys(converter_dict)))])
  l2 = sort(collect(keys(dict_of_filtered_quadratic_elements)))
  @req l1 == l2 "Inconsistency found"

  # (8) We only care about the SR-ideal gens of degree 2. Above, we took care of all relations,
  # (8) for which both variables are not replaced by one of the linear relations. So, let us identify
  # (8) all remaining relations of the SR-ideal, and apply the linear relations to them.
  remaining_relations = Vector{Vector{QQFieldElem}}()
  for my_tuple in ignored_sets

    # The generator must have degree 2 and at least one variable is to be replaced
    if length(my_tuple) == 2 && (my_tuple[1] in bad_positions || my_tuple[2] in bad_positions)

      # Represent first variable by list of coefficients, after plugging in the linear relation
      var1 = zeros(QQ, ncols(my_mat))
      var1[my_tuple[1]] = 1
      if my_tuple[1] in bad_positions
        var1 = lin_rels[my_tuple[1]]
      end

      # Represent second variable by list of coefficients, after plugging in the linear relation
      var2 = zeros(QQ, ncols(my_mat))
      var2[my_tuple[2]] = 1
      if my_tuple[2] in bad_positions
        var2 = lin_rels[my_tuple[2]]
      end

      # Compute the product of the two variables, which represents the new relation
      prod = zeros(QQ, N_filtered_quadratic_elements)
      for k in 1:length(var1)
        if var1[k] != 0
          for l in 1:length(var2)
            if var2[l] != 0
              my_tuple = (min(k, l), max(k, l))
              if haskey(dict_of_filtered_quadratic_elements, my_tuple)
                prod[dict_of_filtered_quadratic_elements[my_tuple]] += var1[k] * var2[l]
              end
            end
          end
        end
      end

      # Remember the result
      push!(remaining_relations, prod)

    end

  end

  # (9) Identify variables that we can remove with the remaining relations
  new_good_positions = 1:N_filtered_quadratic_elements
  if length(remaining_relations) != 0
    remaining_relations_matrix = matrix(QQ, remaining_relations)
    r, new_mat = rref(remaining_relations_matrix)
    @req r == nrows(remaining_relations_matrix) "Cannot remove a variable via linear relations - weird!"
    new_bad_positions = [findfirst(!iszero, row) for row in eachrow(new_mat)]
    new_good_positions = setdiff(1:N_filtered_quadratic_elements, new_bad_positions)
  end

  # (10) Some of the remaining variables are replaced by the final remaining variables
  # (10) Above, we identified the remaining variables. Now we identify how the other variables
  # (10) that appear in dict_of_filtered_quadratic_elements are replaced.
  if length(remaining_relations) != 0
    new_basis = collect(keys(dict_of_filtered_quadratic_elements))
    for (key, value) in dict_of_filtered_quadratic_elements
      if (value in new_good_positions) == false

        # Find relation to repalce this basis element by
        tuple_to_be_replaced = key
        index_of_element_to_be_replaced = value
        row_that_defines_relation = findfirst(k -> k == 1, new_mat[:,index_of_element_to_be_replaced])
        applicable_relation = new_mat[row_that_defines_relation, :]
        applicable_relation[index_of_element_to_be_replaced] = 0
        relation_to_be_applied = (-1) * applicable_relation
        relation_to_be_applied = [[relation_to_be_applied[ivalue], ikey] for (ikey, ivalue) in dict_of_filtered_quadratic_elements if relation_to_be_applied[ivalue] != 0]
        if length(relation_to_be_applied) == 0
          relation_to_be_applied = 0
        end

        # Apply this relation throughout converter_dict, so that this tuple is never used in the values
        for (ikey, ivalue) in converter_dict

          # If the entry maps to zero or is not yet specified, then nothing is to be done. Continue!
          if ivalue == 0 || ivalue === nothing
            continue
          end

          # Is replacement needed?
          tuple_list = [k[2] for k in ivalue]
          position_of_key = findfirst(k -> k == key, tuple_list)
          if position_of_key !== nothing

            # Prepare new lists for the tuples and coefficients
            new_tuple_list = copy(tuple_list)
            new_coeff_list = [k[1] for k in ivalue]

            # Extract the coefficient of interest
            coeff_in_question = new_coeff_list[position_of_key]

            # Remove the tuple to be replaced, and its corresponding coefficient
            deleteat!(new_tuple_list, position_of_key)
            deleteat!(new_coeff_list, position_of_key)

            # Is the list empty after the removal? If so, we map to zero
            if length(new_tuple_list) == 0 && length(new_coeff_list) == 0
              converter_dict[ikey] = 0
              continue
            end

            # Apply relation for element in question
            if relation_to_be_applied == 0
              converter_dict[ikey] = [[new_coeff_list[a], new_tuple_list[a]] for a in 1:length(new_tuple_list)]
            else
              for a in 1:length(relation_to_be_applied)
                if relation_to_be_applied[a][2] in new_tuple_list
                  position_of_tuple = findfirst(k -> k == relation_to_be_applied[a][2], new_tuple_list)
                  new_coeff_list[position_of_tuple] += relation_to_be_applied[a][1]
                else
                  push!(new_tuple_list, relation_to_be_applied[a][2])
                  push!(new_coeff_list, relation_to_be_applied[a][1])
                end
              end
              converter_dict[ikey] == [[new_coeff_list[a], new_tuple_list[a]] for a in 1:length(new_tuple_list)]
            end

          end
        end

        # Assign this value to the tuple to be replaced
        converter_dict[key] = relation_to_be_applied

      end
    end
  end

  # (11) Return the basis elements in terms of cohomology classes
  S = cohomology_ring(v, check = check)
  c_ds = [k.f for k in gens(S)]
  final_list_of_tuples = Tuple{Int64, Int64}[]
  for (key, value) in dict_of_filtered_quadratic_elements
    if value in new_good_positions
      push!(final_list_of_tuples, key)
    end
  end
  basis_of_h22 = [cohomology_class(v, MPolyQuoRingElem(c_ds[my_tuple[1]]*c_ds[my_tuple[2]], S)) for my_tuple in final_list_of_tuples]
  set_attribute!(v, :basis_of_h22_indices, final_list_of_tuples)

  # (12) Consistency check
  l1 = sort(collect(keys(converter_dict))[findall(k -> converter_dict[k] === nothing, collect(keys(converter_dict)))])
  @req l1 == sort(final_list_of_tuples) "Inconsistency found"

  # (13) Convert all entries in converter_dict to cohomology classes and save as attribute
  final_converter_dict = Dict{Tuple{Int, Int}, CohomologyClass}()
  for (key, value) in converter_dict
    if value == nothing
      final_converter_dict[key] = CohomologyClass(v, cohomology_ring(v)(crg[key[1]] *crg[key[2]]))
    elseif value == 0
      final_converter_dict[key] = CohomologyClass(v, zero(cohomology_ring(v, check = check)))
    else
      poly = sum(t[1] * crg[t[2][1]] * crg[t[2][2]] for t in value)
      image_as_cohomology_class = CohomologyClass(v, cohomology_ring(v)(poly))
      final_converter_dict[key] = image_as_cohomology_class
    end
  end
  set_attribute!(v, :converter_dict_h22, final_converter_dict)

  # (14) Return the result - finally!
  return basis_of_h22
end


# The following is an internal function, that is being used to identify all well-quantized G4-fluxes.
# Let G4 in H^(2,2)(toric_ambient_space) a G4-flux ambient space candidate, i.e. the physically
# truly relevant quantity is the restriction of G4 to a hypersurface V(pt) in the toric_ambient_space.
# To tell if this candidate is well-quantized, we execute a necessary check, namely we verify that
# the integral of G4 * [pt] * [d1] * [d2] over X_Sigma is an integer for any two toric divisors d1, d2.
# While we wish to execute this test for any two toric divisors d1, d2 some pairs can be ignored.
# Say, because their intersection locus is trivial because of the SR-ideal, or because their intersection
# has empty intersection with the hypersurface. The following method identifies the remaining pairs of
# toric divisors d1, d2 that we must consider.

@attr Vector{Tuple{Int64, Int64}} function _ambient_space_divisor_pairs_to_be_considered(m::AbstractFTheoryModel)
  gS = gens(cox_ring(ambient_space(m)))
  mnf = Oscar._minimal_nonfaces(ambient_space(m))
  ignored_sets = Set([Tuple(sort(Vector{Int}(Polymake.row(mnf, i)))) for i in 1:Polymake.nrows(mnf)])

  list_of_elements = Vector{Tuple{Int64, Int64}}()
  for k in 1:n_rays(ambient_space(m))
    for l in k:n_rays(ambient_space(m))

      # V(x_k, x_l) = emptyset?
      (k,l) in ignored_sets && continue

      # Simplify the hypersurface polynomial by setting relevant variables to zero.
      # If all coefficients of this new polynomial add to sum, then we keep this generator.
      new_p_hyper = divrem(hypersurface_equation(m), gS[k])[2]
      if k != l
        new_p_hyper = divrem(new_p_hyper, gS[l])[2]
      end
      if sum(coefficients(new_p_hyper)) == 0
        push!(list_of_elements, (k,l))
        continue
      end
      
      # Determine remaining variables, after scaling "away" others.
      remaining_vars_list = Set(1:length(gS))
      for my_exps in ignored_sets
        len_my_exps = length(my_exps)
        inter_len = count(idx -> idx in [k,l], my_exps)
        if (len_my_exps == 2 && inter_len == 1) || (len_my_exps == 3 && inter_len == 2)
          delete!(remaining_vars_list, my_exps[findfirst(idx -> !(idx in [k,l]), my_exps)])
        end
      end
      remaining_vars_list = collect(remaining_vars_list)
  
      # If one monomial of `new_p_hyper` has unset positions (a.k.a. new_p_hyper is not constant upon
      # scaling the remaining variables), then keep this generator.
      for exps in exponents(new_p_hyper)
        if any(x -> x != 0, exps[remaining_vars_list])
          push!(list_of_elements, (k,l))
          break
        end
      end

    end
  end
  
  # Return the findings.
  return list_of_elements
end


# The following is an internal function, that is being used to identify all well-quantized and
# vertical G4-flux ambient candidates. For this, we look for pairs of (pushforwards of) base divisors,
# s.t. their common zero locus does not intersect the CY hypersurface trivially.
# This method makes a pre-selection of such base divisor pairs. "Pre" means that we execute a sufficient,
# but not necessary, check to tell if a pair of base divisors restricts trivially.

@attr Vector{Tuple{Int64, Int64}} function _ambient_space_base_divisor_pairs_to_be_considered(m::AbstractFTheoryModel)
  gS = gens(cox_ring(ambient_space(m)))
  mnf = Oscar._minimal_nonfaces(ambient_space(m))
  ignored_sets = Set([Tuple(sort(Vector{Int}(Polymake.row(mnf, i)))) for i in 1:Polymake.nrows(mnf)])

  list_of_elements = Vector{Tuple{Int64, Int64}}()
  for k in 1:n_rays(base_space(m))
    for l in k:n_rays(base_space(m))

      # V(x_k, x_l) = emptyset?
      (k,l) in ignored_sets && continue

      # Simplify the hypersurface polynomial by setting relevant variables to zero.
      # If all coefficients of this new polynomial add to sum, then we keep this generator.
      new_p_hyper = divrem(hypersurface_equation(m), gS[k])[2]
      if k != l
        new_p_hyper = divrem(new_p_hyper, gS[l])[2]
      end
      if sum(coefficients(new_p_hyper)) == 0
        push!(list_of_elements, (k,l))
        continue
      end
      
      # Determine remaining variables, after scaling "away" others.
      remaining_vars_list = Set(1:length(gS))
      for my_exps in ignored_sets
        len_my_exps = length(my_exps)
        inter_len = count(idx -> idx in [k,l], my_exps)
        if (len_my_exps == 2 && inter_len == 1) || (len_my_exps == 3 && inter_len == 2)
          delete!(remaining_vars_list, my_exps[findfirst(idx -> !(idx in [k,l]), my_exps)])
        end
      end
      remaining_vars_list = collect(remaining_vars_list)
  
      # If one monomial of `new_p_hyper` has unset positions (a.k.a. new_p_hyper is not constant upon
      # scaling the remaining variables), then keep this generator.
      for exps in exponents(new_p_hyper)
        if any(x -> x != 0, exps[remaining_vars_list])
          push!(list_of_elements, (k,l))
          break
        end
      end

    end
  end
  
  # Return the findings.
  return list_of_elements
end


function _ambient_space_models_of_g4_fluxes(m::AbstractFTheoryModel; check::Bool = true)

  # Entry check
  @req base_space(m) isa NormalToricVariety "Base space must be a toric variety for computation of ambient space G4 candidates"
  if has_attribute(m, :ambient_space_models_of_g4_fluxes)
    return get_attribute(m, :ambient_space_models_of_g4_fluxes)::Vector{CohomologyClass}
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
  return filtered_h22_basis::Vector{CohomologyClass}

end


@doc raw"""
    chosen_g4_flux_basis(m::AbstractFTheoryModel; check::Bool = true)::Vector{CohomologyClass}

Given an F-theory model `m` defined as a hypersurface in a simplicial and
complete toric base, this method computes a basis of $H^{2,2}(X, \mathbb{Q})$
(using the method `basis_of_h22`) and then filters out certain basis elements 
whose restriction to the hypersurface in question is trivial. The criteria for 
"certain" elements are explained in the documentation above this method.

Note: Checking whether a toric variety $X$ is complete and simplicial can be
computationally expensive. The optional argument `check` can be set to `false`
to skip these tests.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> Kbar = anticanonical_divisor_class(B3)
Divisor class on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w"=>Kbar))
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> g4_basis = chosen_g4_flux_basis(t);

julia> length(g4_basis)
2

julia> g4_basis[1]
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not executed

julia> cohomology_class(g4_basis[1])
Cohomology class on a normal toric variety given by z^2

julia> cohomology_class(g4_basis[2])
Cohomology class on a normal toric variety given by y^2

julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 8))
Hypersurface model over a concrete base

julia> g4_basis = chosen_g4_flux_basis(qsm_model, check = false);

julia> cohomology_class(g4_basis[1])
Cohomology class on a normal toric variety given by x15*e2

julia> length(g4_basis) == 172
true
```
"""
chosen_g4_flux_basis(m::AbstractFTheoryModel; check::Bool = true) = [G4Flux(m, c) for c in _ambient_space_models_of_g4_fluxes(m, check = check)]
