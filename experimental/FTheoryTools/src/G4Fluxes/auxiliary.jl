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
  set_attribute!(v, :basis_of_h22, basis_of_h22)
  set_attribute!(v, :basis_of_h22_indices, final_list_of_tuples)
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

function _ambient_space_divisor_pairs_to_be_considered(m::AbstractFTheoryModel)::Vector{Tuple{Int64, Int64}}

  if has_attribute(m, :_ambient_space_divisor_pairs_to_be_considered)
    return get_attribute(m, :_ambient_space_divisor_pairs_to_be_considered)
  end

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
  
  # Remember this result as attribute and return the findings.
  set_attribute!(m, :_ambient_space_divisor_pairs_to_be_considered, list_of_elements)
  return list_of_elements
  
end


# The following is an internal function, that is being used to identify all well-quantized and
# vertical G4-flux ambient candidates. For this, we look for pairs of (pushforwards of) base divisors,
# s.t. their common zero locus does not intersect the CY hypersurface trivially.
# This method makes a pre-selection of such base divisor pairs. "Pre" means that we execute a sufficient,
# but not necessary, check to tell if a pair of base divisors restricts trivially.

function _ambient_space_base_divisor_pairs_to_be_considered(m::AbstractFTheoryModel)::Vector{Tuple{Int64, Int64}}

  if has_attribute(m, :_ambient_space_base_divisor_pairs_to_be_considered)
    return get_attribute(m, :_ambient_space_base_divisor_pairs_to_be_considered)
  end

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
  
  # Remember this result as attribute and return the findings.
  set_attribute!(m, :_ambient_space_base_divisor_pairs_to_be_considered, list_of_elements)
  return list_of_elements
  
end
