###########################################################################################################################
# The following code deals with finding a basis of H22 of the ambient toric space, respectively converting into this basis.
# The following code deals with finding a basis of H22 of the ambient toric space, respectively converting into this basis.
###########################################################################################################################


@doc raw"""
    converter_dict_h22_ambient(m::AbstractFTheoryModel; check::Bool = true)

By virtue of Theorem 12.4.1 in [CLS11](@cite), one can compute a monomial
basis of ``H^4(X, \mathbb{Q})`` for a simplicial, complete toric variety ``X``
by truncating its cohomology ring to degree ``2``. Inspired by this, this
method identifies a basis of ``H^{(2,2)}(X, \mathbb{Q})`` by multiplying
pairs of cohomology classes associated with toric coordinates.

By definition, ``H^{(2,2)}(X, \mathbb{Q})`` is a subset of ``H^{4}(X, \mathbb{Q})``.
However, by Theorem 9.3.2 in [CLS11](@cite), for complete and simplicial
toric varieties and ``p \neq q`` it holds ``H^{(p,q)}(X, \mathbb{Q}) = 0``. It follows
that for such varieties ``H^{(2,2)}(X, \mathbb{Q}) = H^4(X, \mathbb{Q})`` and the
vector space dimension of those spaces agrees with the Betti number ``b_4(X)``.

This method computes a dictionary which allows to map any element of
``H^(2,2)(X, \mathbb{Q})`` into the basis found/chosen for ``H^(2,2)(X, \mathbb{Q})``.

Note that it can be computationally very demanding to check if a toric variety
``X`` is complete (and simplicial). The optional argument `check` can be set
to `false` to skip these tests.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 283))
Hypersurface model over a concrete base

julia> cdh22 = converter_dict_h22_ambient(qsm_model, check = false);

julia> length(collect(keys(cdh22)))
81
```
"""
@attr Dict{Tuple{Int64, Int64}, Vector{Tuple{QQFieldElem, Tuple{Int64, Int64}}}} function converter_dict_h22_ambient(m::AbstractFTheoryModel; check::Bool = true)

  # (1) Entry checks
  @req base_space(m) isa NormalToricVariety "Computation of converter_dict_h22_ambient only supported for toric base and ambient spaces"
  @req dim(ambient_space(m)) == 5 "Computation of converter_dict_h22_ambient only supported for 5-dimensional toric ambient spaces"
  if check
    @req is_complete(ambient_space(m)) "Computation of converter_dict_h22_ambient only supported for complete toric ambient spaces"
    @req is_simplicial(ambient_space(m)) "Computation of converter_dict_h22_ambient only supported for simplicial toric ambient space"
  end


  # (2) Initialize the dict that converts the "naive" generating set into our chosen basis.
  # (2) Map pairs, that are zero by virtue of the Stanley-Reisner ideal to zero.
  v = ambient_space(m)
  mnf = Oscar._minimal_nonfaces(v)
  ignored_sets = Set([Tuple(sort(Vector{Int}(Polymake.row(mnf, i)))) for i in 1:Polymake.nrows(mnf)])
  converter_dict = Dict{Tuple{Int, Int}, Vector{Tuple{QQFieldElem, Tuple{Int64, Int64}}}}()
  for k in 1:n_rays(v)
    for l in k:n_rays(v)
      my_tuple = (k, l)
      if (my_tuple in ignored_sets)
        converter_dict[my_tuple] = Vector{Tuple{QQFieldElem, Tuple{Int64, Int64}}}()
      else
        converter_dict[my_tuple] = [(QQ(1), my_tuple)]
      end
    end
  end


  # (3) With the linear relations, we can remove some variables. Identify these relations.
  # (3) Also, find the remaining "good" variables: good_positions encodes the variables that we keep. Note that this is a basis of H^(1,1).
  my_mat = matrix(QQ, transpose(matrix(ZZ, rays(v))))
  N_lin_rel, my_mat = rref(my_mat)
  @req N_lin_rel == nrows(my_mat) "Inconsistency detected - cannot remove as many variables as there are linear relations"
  bad_positions = [findfirst(!iszero, row) for row in eachrow(my_mat)]
  lin_rels = Dict{Int, Vector{QQFieldElem}}()
  for k in 1:nrows(my_mat)
    my_relation = (-1) * my_mat[k, :]
    my_relation[bad_positions[k]] = 0
    @req all(k -> k == 0, my_relation[bad_positions]) "Inconsistency!"
    lin_rels[bad_positions[k]] = my_relation
  end
  good_positions = setdiff(1:n_rays(v), bad_positions)
  

  # (4) Apply linear relations to converter_dict, to eliminate the variables in bad_positions.
  for (key, value) in converter_dict
    if value != 0
      image = Vector{Tuple{QQFieldElem, Tuple{Int64, Int64}}}()

      # Only the first entry needs replacing by the linear relations
      if (key[1] in bad_positions) && (key[2] in bad_positions) == false
        my_relation = lin_rels[key[1]]
        posi = findall(k -> k != 0, my_relation)
        coeffs = my_relation[posi]
        for k in 1:length(posi)
          new_tuple = (min(posi[k], key[2]), max(posi[k], key[2]))
          if !(new_tuple in ignored_sets)
            current_tuples = [k[2] for k in image]
            position_of_interest = findfirst(d -> d == new_tuple, current_tuples)
            if position_of_interest === nothing
              push!(image, (coeffs[k], new_tuple))
            else
              old_entry = image[position_of_interest]
              new_entry = (old_entry[1] + coeffs[k], old_entry[2])
              image[position_of_interest] = new_entry
            end
          end
        end
      end
      
      # Only the second entry needs replacing by the linear relations
      if (key[2] in bad_positions) && (key[1] in bad_positions) == false
        my_relation = lin_rels[key[2]]
        posi = findall(k -> k != 0, my_relation)
        coeffs = my_relation[posi]
        for k in 1:length(posi)
          new_tuple = (min(key[1], posi[k]), max(key[1], posi[k]))
          if !(new_tuple in ignored_sets)
            current_tuples = [k[2] for k in image]
            position_of_interest = findfirst(d -> d == new_tuple, current_tuples)
            if position_of_interest === nothing
              push!(image, (coeffs[k], new_tuple))
            else
              old_entry = image[position_of_interest]
              new_entry = (old_entry[1] + coeffs[k], old_entry[2])
              image[position_of_interest] = new_entry
            end
          end
        end
      end

      # Both entries need replacing by the linear relations
      if (key[1] in bad_positions) && (key[2] in bad_positions)
        my_relation = lin_rels[key[1]]
        posi = findall(k -> k != 0, my_relation)
        coeffs = my_relation[posi]
        for k in 1:length(posi)
          new_tuple = (posi[k], key[2])
          current_tuples = [k[2] for k in image]
          position_of_interest = findfirst(d -> d == new_tuple, current_tuples)
          if position_of_interest === nothing
            push!(image, (coeffs[k], new_tuple))
          else
            image[position_of_interest][1] += coeffs[k]
          end
        end
        my_relation = lin_rels[key[2]]
        posi = findall(k -> k != 0, my_relation)
        coeffs = my_relation[posi]
        old_image = copy(image)
        image = Vector{Tuple{QQFieldElem, Tuple{Int64, Int64}}}()
        for i in 1:length(old_image)
          for k in 1:length(posi)
            old_coeff = old_image[i][1]
            new_tuple = (min(old_image[i][2][1], posi[k]), max(old_image[i][2][1], posi[k]))
            if !(new_tuple in ignored_sets)
              current_tuples = [k[2] for k in image]
              position_of_interest = findfirst(d -> d == new_tuple, current_tuples)
              if position_of_interest === nothing
                push!(image, (old_coeff * coeffs[k], new_tuple))
              else
                old_entry = image[position_of_interest]
                new_entry = (old_entry[1] + old_coeff * coeffs[k], old_entry[2])
                image[position_of_interest] = new_entry
              end
            end
          end
        end
      end
      
      # If there was a replacement, update the value in the dict
      if length(image) > 0
        converter_dict[key] = image
      end
    end
  end


  # (5) There can be left-over relations after applying the above. We proceed as follows:
  # (5) Step 1: Make a list of all tuples (a,b) of divisor pairs that need considering after the previous step. We order them, i.e. assign a number to each tuple (a,b).
  # (5) Step 2: Identify the remaining relations among these elements.
  # (5) Step 3: Identify which of the tuples determined in step 1 can be replaced by Q-linear relations of the others.

  # (5a) Step 1
  N_filtered_quadratic_elements = 0
  dict_of_filtered_quadratic_elements = Dict{Tuple{Int64, Int64}, Int64}()
  for k in 1:length(good_positions)
    for l in k:length(good_positions)
      my_tuple = (min(good_positions[k], good_positions[l]), max(good_positions[k], good_positions[l]))
      if !(my_tuple in ignored_sets)
        N_filtered_quadratic_elements += 1
        dict_of_filtered_quadratic_elements[my_tuple] = N_filtered_quadratic_elements
      end
    end
  end
  # Consistency check
  a1 = sort(unique(vcat([[u[2] for u in k] for k in filter(x -> x != 0, collect(values(converter_dict)))]...)))
  a2 = sort(collect(keys(dict_of_filtered_quadratic_elements)))
  @req a1 == a2 "Inconsistency encountered"


  # (5b) Step 2
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


  # (5c) Identify variables that we can remove with the remaining relations, then remove them.
  #new_good_positions = 1:N_filtered_quadratic_elements
  if length(remaining_relations) != 0
    remaining_relations_matrix = matrix(QQ, remaining_relations)
    r, new_mat = rref(remaining_relations_matrix)
    @req r == nrows(remaining_relations_matrix) "Cannot remove a variable via linear relations - weird!"
    new_bad_positions = [findfirst(!iszero, row) for row in eachrow(new_mat)]
    new_good_positions = setdiff(1:N_filtered_quadratic_elements, new_bad_positions)

    # Consistency check
    basis_elements = Tuple{Int, Int}[]
    for k in new_good_positions
      matching_keys = [key for (key, value) in dict_of_filtered_quadratic_elements if value == k]    
      if length(matching_keys) != 1
        println("Verification failed: Expected exactly one key for value $k, but found $(length(matching_keys))")
        return false, []
      end
      push!(basis_elements, matching_keys[1])
    end

    # Apply the remaining relations
    for (key, value) in dict_of_filtered_quadratic_elements
      if (value in new_good_positions) == false

        # Find relation to repalce this basis element by
        row_that_defines_relation = findfirst(k -> k == 1, new_mat[:,value])
        applicable_relation = copy(new_mat[row_that_defines_relation, :])
        applicable_relation[value] = 0
        relation_to_be_applied = (-1) * applicable_relation
        relation_to_be_applied = [(relation_to_be_applied[ivalue], ikey) for (ikey, ivalue) in dict_of_filtered_quadratic_elements if relation_to_be_applied[ivalue] != 0]
        if length(relation_to_be_applied) == 0
          relation_to_be_applied = 0
        end

        # Apply this relation throughout converter_dict, so that this tuple is never used in the values
        for (ikey, ivalue) in converter_dict

          # If the entry maps to zero, then nothing is to be done. Continue!
          if length(ivalue) == 0
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
            if length(new_tuple_list) == 0 && length(new_coeff_list) == 0 && relation_to_be_applied == 0
              converter_dict[ikey] = Vector{Tuple{QQFieldElem, Tuple{Int64, Int64}}}()
            elseif length(new_tuple_list) == 0 && length(new_coeff_list) == 0
              converter_dict[ikey] = [(coeff_in_question * k[1], k[2]) for k in relation_to_be_applied if k[1] != 0]
            elseif relation_to_be_applied == 0
              converter_dict[ikey] = [(new_coeff_list[a], new_tuple_list[a]) for a in 1:length(new_tuple_list)]
            else
              for a in 1:length(relation_to_be_applied)
                if relation_to_be_applied[a][2] in new_tuple_list
                  position_of_tuple = findfirst(k -> k == relation_to_be_applied[a][2], new_tuple_list)
                  new_coeff_list[position_of_tuple] += coeff_in_question * relation_to_be_applied[a][1]
                else
                  push!(new_tuple_list, relation_to_be_applied[a][2])
                  push!(new_coeff_list, coeff_in_question * relation_to_be_applied[a][1])
                end
              end
              converter_dict[ikey] = [(new_coeff_list[a], new_tuple_list[a]) for a in 1:length(new_tuple_list) if new_coeff_list[a] != 0]
            end

          end
        end # End updating converter dict

      end
    end

    # Final consistency check
    a1 = sort(unique(vcat([[u[2] for u in k] for k in filter(x -> x != 0, collect(values(converter_dict)))]...)))
    @req a1 == basis_elements "Inconsistency encountered"
  
  end

  # To save memory, we filter out any tuples k which map to zero
  return Dict(k => val for (k, val) in converter_dict if !isempty(val))

end


@doc raw"""
    basis_of_h22_ambient_indices(m::AbstractFTheoryModel; check::Bool = true)

By virtue of Theorem 12.4.1 in [CLS11](@cite), one can compute a monomial
basis of ``H^4(X, \mathbb{Q})`` for a simplicial, complete toric variety ``X``
by truncating its cohomology ring to degree ``2``. Inspired by this, this
method identifies a basis of ``H^{(2,2)}(X, \mathbb{Q})`` by multiplying
pairs of cohomology classes associated with toric coordinates.

By definition, ``H^{(2,2)}(X, \mathbb{Q})`` is a subset of ``H^{4}(X, \mathbb{Q})``.
However, by Theorem 9.3.2 in [CLS11](@cite), for complete and simplicial
toric varieties and ``p \neq q`` it holds ``H^{(p,q)}(X, \mathbb{Q}) = 0``. It follows
that for such varieties ``H^{(2,2)}(X, \mathbb{Q}) = H^4(X, \mathbb{Q})`` and the
vector space dimension of those spaces agrees with the Betti number ``b_4(X)``.

This method computes a vector of tuples of two integers. The tuple ``(a,b)`` states
that the product of the a-th and b-th variables in the Cox ring is an element of
the basis of ``H^(2,2)(X, \mathbb{Q})``, that we have chosen.

Note that it can be computationally very demanding to check if a toric variety
``X`` is complete (and simplicial). The optional argument `check` can be set
to `false` to skip these tests.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 283))
Hypersurface model over a concrete base

julia> bhi = basis_of_h22_ambient_indices(qsm_model, check = false);

julia> length(bhi) == betti_number(ambient_space(qsm_model), 4)
true
```
"""
@attr Vector{Tuple{Int64, Int64}} function basis_of_h22_ambient_indices(m::AbstractFTheoryModel; check::Bool = true)
  cdh22 = converter_dict_h22_ambient(m, check = check)
  return unique(vcat([[u[2] for u in k] for k in filter(x -> x != 0, collect(values(cdh22)))]...))
end


@doc raw"""
    basis_of_h22_ambient(m::AbstractFTheoryModel; check::Bool = true)

By virtue of Theorem 12.4.1 in [CLS11](@cite), one can compute a monomial
basis of ``H^4(X, \mathbb{Q})`` for a simplicial, complete toric variety ``X``
by truncating its cohomology ring to degree ``2``. Inspired by this, this
method identifies a basis of ``H^{(2,2)}(X, \mathbb{Q})`` by multiplying
pairs of cohomology classes associated with toric coordinates.

By definition, ``H^{(2,2)}(X, \mathbb{Q})`` is a subset of ``H^{4}(X, \mathbb{Q})``.
However, by Theorem 9.3.2 in [CLS11](@cite), for complete and simplicial
toric varieties and ``p \neq q`` it holds ``H^{(p,q)}(X, \mathbb{Q}) = 0``. It follows
that for such varieties ``H^{(2,2)}(X, \mathbb{Q}) = H^4(X, \mathbb{Q})`` and the
vector space dimension of those spaces agrees with the Betti number ``b_4(X)``.

This method computes a basis of ``H^{(2,2)}(X, \mathbb{Q})`` as vector of cohomology
classes.

Note that it can be computationally very demanding to check if a toric variety
``X`` is complete (and simplicial). The optional argument `check` can be set
to `false` to skip these tests.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 283))
Hypersurface model over a concrete base

julia> length(basis_of_h22_ambient(qsm_model, check = false)) == betti_number(ambient_space(qsm_model), 4)
true
```
"""
@attr Vector{CohomologyClass} function basis_of_h22_ambient(m::AbstractFTheoryModel; check::Bool = true)
  basis_indices = basis_of_h22_ambient_indices(m, check = false)
  v = ambient_space(m)
  S = cohomology_ring(v, check = check)
  c_ds = [lift(k) for k in gens(S)]
  return [cohomology_class(v, MPolyQuoRingElem(c_ds[mt[1]]*c_ds[mt[2]], S)) for mt in basis_indices]
end



###########################################################################################################################
# The following code deals with finding a basis of H22 of the hypersurface, respectively converting into this basis.
# The following code deals with finding a basis of H22 of the hypersurface, respectively converting into this basis.
###########################################################################################################################

@doc raw"""
    gens_of_h22_hypersurface_indices(m::AbstractFTheoryModel; check::Bool = true)

By virtue of Theorem 12.4.1 in [CLS11](@cite), one can compute a monomial
basis of ``H^4(X, \mathbb{Q})`` for a simplicial, complete toric variety ``X``
by truncating its cohomology ring to degree ``2``. Inspired by this, this
method identifies a basis of ``H^{(2,2)}(X, \mathbb{Q})`` by multiplying
pairs of cohomology classes associated with toric coordinates.

By definition, ``H^{(2,2)}(X, \mathbb{Q})`` is a subset of ``H^{4}(X, \mathbb{Q})``.
However, by Theorem 9.3.2 in [CLS11](@cite), for complete and simplicial
toric varieties and ``p \neq q`` it holds ``H^{(p,q)}(X, \mathbb{Q}) = 0``. It follows
that for such varieties ``H^{(2,2)}(X, \mathbb{Q}) = H^4(X, \mathbb{Q})`` and the
vector space dimension of those spaces agrees with the Betti number ``b_4(X)``.

Once restricted to the hypersurface of the model, these elements are generators of
a subset ``S \subseteq H^(2,2)(Y, \mathbb{Q})``, where ``Y`` is the hypersurface of
interest in ``X``. This method computes a vector of tuples, each tuple encodes on
generator of ``S \subseteq H^(2,2)(Y, \mathbb{Q})``. Specifically, the tuple ``(a,b)``
indicates that the product of the ``a``-th and the ``b``-th variable in the Cox ring
encodes an element of the cohomology ring of ``X``, whose restriction to ``Y`` is one
of our generators.

Note that it can be computationally very demanding to check if a toric variety
``X`` is complete (and simplicial). The optional argument `check` can be set
to `false` to skip these tests.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 283))
Hypersurface model over a concrete base

julia> length(gens_of_h22_hypersurface_indices(qsm_model, check = false))
25
```
"""
@attr Vector{Tuple{Int64, Int64}} function gens_of_h22_hypersurface_indices(m::AbstractFTheoryModel; check::Bool = true)

  filtered_h22_basis_indices = deepcopy(basis_of_h22_ambient_indices(m, check = check))
  gS = gens(cox_ring(ambient_space(m)))
  mnf = Oscar._minimal_nonfaces(ambient_space(m))
  sr_ideal_pos = [Vector{Int}(Polymake.row(mnf, i)) for i in 1:Polymake.nrows(mnf)]
  for a in length(filtered_h22_basis_indices):-1:1
    
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
      deleteat!(filtered_h22_basis_indices, a)
    end

  end
  return filtered_h22_basis_indices

end

@doc raw"""
    gens_of_h22_hypersurface(m::AbstractFTheoryModel; check::Bool = true)

By virtue of Theorem 12.4.1 in [CLS11](@cite), one can compute a monomial
basis of ``H^4(X, \mathbb{Q})`` for a simplicial, complete toric variety ``X``
by truncating its cohomology ring to degree ``2``. Inspired by this, this
method identifies a basis of ``H^{(2,2)}(X, \mathbb{Q})`` by multiplying
pairs of cohomology classes associated with toric coordinates.

By definition, ``H^{(2,2)}(X, \mathbb{Q})`` is a subset of ``H^{4}(X, \mathbb{Q})``.
However, by Theorem 9.3.2 in [CLS11](@cite), for complete and simplicial
toric varieties and ``p \neq q`` it holds ``H^{(p,q)}(X, \mathbb{Q}) = 0``. It follows
that for such varieties ``H^{(2,2)}(X, \mathbb{Q}) = H^4(X, \mathbb{Q})`` and the
vector space dimension of those spaces agrees with the Betti number ``b_4(X)``.

Once restricted to the hypersurface of the model, these elements are generators of
a subset ``S \subseteq H^(2,2)(Y, \mathbb{Q})``, where ``Y`` is the hypersurface of
interest in ``X``. This method computes this generating set in terms of ambient
space cohomology classes.

Note that it can be computationally very demanding to check if a toric variety
``X`` is complete (and simplicial). The optional argument `check` can be set
to `false` to skip these tests.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 283))
Hypersurface model over a concrete base

julia> length(gens_of_h22_hypersurface(qsm_model, check = false))
25
```
"""
@attr Vector{CohomologyClass} function gens_of_h22_hypersurface(m::AbstractFTheoryModel; check::Bool = true)
  basis_indices = gens_of_h22_hypersurface_indices(m, check = false)
  v = ambient_space(m)
  S = cohomology_ring(v, check = check)
  c_ds = [lift(k) for k in gens(S)]
  return [cohomology_class(v, MPolyQuoRingElem(c_ds[mt[1]]*c_ds[mt[2]], S)) for mt in basis_indices]
end


@doc raw"""
    converter_dict_h22_hypersurface(m::AbstractFTheoryModel; check::Bool = true)

By virtue of Theorem 12.4.1 in [CLS11](@cite), one can compute a monomial
basis of ``H^4(X, \mathbb{Q})`` for a simplicial, complete toric variety ``X``
by truncating its cohomology ring to degree ``2``. Inspired by this, this
method identifies a basis of ``H^{(2,2)}(X, \mathbb{Q})`` by multiplying
pairs of cohomology classes associated with toric coordinates.

By definition, ``H^{(2,2)}(X, \mathbb{Q})`` is a subset of ``H^{4}(X, \mathbb{Q})``.
However, by Theorem 9.3.2 in [CLS11](@cite), for complete and simplicial
toric varieties and ``p \neq q`` it holds ``H^{(p,q)}(X, \mathbb{Q}) = 0``. It follows
that for such varieties ``H^{(2,2)}(X, \mathbb{Q}) = H^4(X, \mathbb{Q})`` and the
vector space dimension of those spaces agrees with the Betti number ``b_4(X)``. Once
restricted to the hypersurface ``Y`` in question, we obtain a generating set of a subset
``S \subseteq H^{(2,2)}(Y, \mathbb{Q})``.

This method computes a dictionary which allows to map any element of ``H^(2,2)(X, \mathbb{Q})``
into our chosen generating set of ``S \subseteq H^{(2,2)}(Y, \mathbb{Q})``.

Note that it can be computationally very demanding to check if a toric variety
``X`` is complete (and simplicial). The optional argument `check` can be set
to `false` to skip these tests.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)))
julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 283))
Hypersurface model over a concrete base

julia> cdh22 = converter_dict_h22_hypersurface(qsm_model, check = false);

julia> length(collect(keys(cdh22)))
81
```
"""
@attr Dict{Tuple{Int64, Int64}, Vector{Tuple{QQFieldElem, Tuple{Int64, Int64}}}} function converter_dict_h22_hypersurface(m::AbstractFTheoryModel; check::Bool = true)
  non_trivial_indices = gens_of_h22_hypersurface_indices(m, check = check)
  old_indices = basis_of_h22_ambient_indices(m, check = check)
  to_be_deleted = setdiff(old_indices, non_trivial_indices)
  new_converter = deepcopy(converter_dict_h22_ambient(m, check = check))
  for k in 1:length(to_be_deleted)
    for (key, value) in new_converter
      tuple_list = [a[2] for a in value]
      if to_be_deleted[k] in tuple_list
        # There is a tuple to be deleted...
        posi = findfirst(a -> a == to_be_deleted[k], tuple_list)
        new_value = deleteat!(value, posi)
        new_converter[key] = new_value
      end
    end
  end
  # To save memory, we filter out any tuples k which map to zero
  return Dict(k => v for (k, v) in new_converter if !isempty(v))
end


@doc raw"""
    chosen_g4_flux_gens(m::AbstractFTheoryModel; check::Bool = true)::Vector{CohomologyClass}

Given an F-theory model `m` defined as a hypersurface in a simplicial and
complete toric base, this method computes a basis of ``H^{2,2}(X, \mathbb{Q})``
(using the method `basis_of_h22`) and then filters out certain basis elements 
whose restriction to the hypersurface in question is trivial. The criteria for 
"certain" elements are explained in the documentation above this method.

Note: Checking whether a toric variety ``X`` is complete and simplicial can be
computationally expensive. The optional argument `check` can be set to `false`
to skip these tests.

# Examples
```jldoctest; setup = :(Oscar.LazyArtifacts.ensure_artifact_installed("QSMDB", Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir))), filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> Kbar = anticanonical_divisor_class(B3)
Divisor class on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w"=>Kbar))
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> g4_basis = chosen_g4_flux_gens(t);

julia> length(g4_basis)
2

julia> g4_basis[1]
G4-flux candidate
  - Elementary quantization checks: not executed
  - Transversality checks: not executed
  - Non-abelian gauge group: breaking pattern not analyzed
  - Tadpole cancellation check: not computed

julia> cohomology_class(g4_basis[1])
Cohomology class on a normal toric variety given by y^2

julia> cohomology_class(g4_basis[2])
Cohomology class on a normal toric variety given by z^2

julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 8))
Hypersurface model over a concrete base

julia> g4_basis = chosen_g4_flux_gens(qsm_model, check = false);

julia> cohomology_class(g4_basis[1])
Cohomology class on a normal toric variety given by x15*e2

julia> length(g4_basis) == 172
true
```
"""
@attr Vector{G4Flux} function chosen_g4_flux_gens(m::AbstractFTheoryModel; check::Bool = true)
  gens = [G4Flux(m, c) for c in gens_of_h22_hypersurface(m, check = check)]
  for k in 1:length(gens)
    set_attribute!(gens[k], :offset, zeros(Int, length(gens)))
    flux_coords = zeros(Int, length(gens))
    flux_coords[k] = 1
    set_attribute!(gens[k], :flux_coordinates, flux_coords)
  end
  return gens
end



###########################################################################################################################
# The following are internal helper functions to find fluxes with specific properties.
# The following are internal helper functions to find fluxes with specific properties.
###########################################################################################################################

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
