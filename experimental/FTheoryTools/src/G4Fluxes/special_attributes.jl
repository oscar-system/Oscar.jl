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
function basis_of_h22(v::NormalToricVariety; check::Bool = true)::Vector{CohomologyClass}

  # (0) Some initial checks
  if check
    @req is_complete(v) "Computation of basis of H22 is currently only supported for complete toric varieties"
    @req is_simplicial(v) "Computation of basis of H22 is currently only supported for simplicial toric varieties"
  end
  if dim(v) < 4
    set_attribute!(v, :basis_of_h22, Vector{CohomologyClass}())
  end
  if has_attribute(v, :basis_of_h22)
    return get_attribute(v, :basis_of_h22)
  end

  # (1) Prepare some data of the variety
  mnf = Oscar._minimal_nonfaces(v)
  ignored_sets = Set([Tuple(sort(Vector{Int}(Polymake.row(mnf, i)))) for i in 1:Polymake.nrows(mnf)])

  # (2) Prepare the linear relations
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

  # (3) Prepare a list of those variables that we keep, a.k.a. a basis of H^(1,1)
  good_positions = setdiff(1:n_rays(v), bad_positions)
  n_good_positions = length(good_positions)

  # (4) Make a list of all quadratic elements in the cohomology ring, which are not generators of the SR-ideal.
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

  # (5) We only care about the SR-ideal gens of degree 2. Above, we took care of all relations,
  # (5) for which both variables are not replaced by one of the linear relations. So, let us identify
  # (5) all remaining relations of the SR-ideal, and apply the linear relations to them.
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
  
  # (10) Return the basis elements in terms of cohomology classes
  S = cohomology_ring(v, check = check)
  c_ds = [k.f for k in gens(S)]
  final_list_of_tuples = []
  for (key, value) in dict_of_filtered_quadratic_elements
    if value in new_good_positions
      push!(final_list_of_tuples, key)
    end
  end
  basis_of_h22 = [cohomology_class(v, MPolyQuoRingElem(c_ds[my_tuple[1]]*c_ds[my_tuple[2]], S)) for my_tuple in final_list_of_tuples]
  #set_attribute!(v, :basis_of_h22, basis_of_h22)
  return basis_of_h22

end


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
  
  # Prepare data of the toric ambient space
  gS = gens(cox_ring(ambient_space(m)))
  mnf = Oscar._minimal_nonfaces(ambient_space(m))
  sr_ideal_pos = [Vector{Int}(Polymake.row(mnf, i)) for i in 1:Polymake.nrows(mnf)]
  filtered_h22_basis = basis_of_h22(ambient_space(m), check = check)

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

  return filtered_h22_basis
end
