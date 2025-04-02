# ---------------------------------------------------------------------------------------------------------
# (1) Compute the intersection product of an algebraic cycle with a hypersurface.
# ---------------------------------------------------------------------------------------------------------

function sophisticated_intersection_product(v::NormalToricVariety, indices::NTuple{4, Int64}, hypersurface_equation::MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}, inter_dict::Dict{NTuple{4, Int64}, ZZRingElem}, s_inter_dict::Dict{String, ZZRingElem}, data::NamedTuple)

  # (A) Have we computed this intersection number in the past? If so, just use that result...
  indices = Tuple(sort(collect(indices)))
  if haskey(inter_dict, indices)
    return inter_dict[indices]
  end
  
  # (B) Is the intersection, by virtue of the SR-ideal, trivial?
  for sr_set in data.sr_ideal_pos
    if is_subset(sr_set, indices)
      inter_dict[indices] = ZZ(0)
      return ZZ(0)
    end
  end

  # (C) Deal with self-intersection and should-never-happen case.
  distinct_variables = Set(indices)
  if length(distinct_variables) < 4 && length(distinct_variables) >= 1
    return intersection_from_equivalent_cycle(v, indices, hypersurface_equation, inter_dict, s_inter_dict, data)
  end
  if length(distinct_variables) == 0
    println("WEIRD! THIS SHOULD NEVER HAPPEN! INFORM THE AUTHORS!")
    println("")
  end


  # (D) Deal with transverse intersection...

  # D.1 Work out the intersection locus in detail.
  pt_reduced, gs_reduced, remaining_vars, reduced_scaling_relations = Oscar._reduce_hypersurface_equation(v, hypersurface_equation, indices, data)

  # D.2 If pt == 0, then we are not looking at a transverse intersection. So take an equivalent cycle and try again...
  if is_zero(pt_reduced)
    return intersection_from_equivalent_cycle(v, indices, hypersurface_equation, inter_dict, s_inter_dict, data)
  end

  # D.3 If pt is constant and non-zero, then the intersection is trivial.
  if is_constant(pt_reduced) && is_zero(pt_reduced) == false
    inter_dict[indices] = ZZ(0)
    return ZZ(0)
  end

  # D.4 Helper function for the cases below
  function has_one_and_rest_zero(vec)
    return count(==(1), vec) == 1 && all(x -> x == 0 || x == 1, vec)
  end

  # C.5 Cover a case that seems to appear frequently for our investigation:
  # pt_reduced of the form a * x + b * y for non-zero number a,b and remaining variables x, y subject to a reduced SR generator x * y and scaling relation [1,1].
  # This will thus always give exactly one solution (x = 1, y = -a/b), and so the intersection number is one.
  if length(gs_reduced) == 1 && length(remaining_vars) == 2
    mons_list = collect(monomials(pt_reduced))
    if length(mons_list) == 2
      if all(x -> x != 0, collect(coefficients(pt_reduced)))
        exps_list = [collect(exponents(k))[1] for k in mons_list]
        if has_one_and_rest_zero(exps_list[1]) && has_one_and_rest_zero(exps_list[2])
          if gs_reduced[1] == remaining_vars[1] * remaining_vars[2]
            if reduced_scaling_relations == matrix(ZZ, [[1,1]])
              inter_dict[indices] = ZZ(1)
              return ZZ(1)
            end
          end
        end
      end
    end
  end

  # C.6 Cover a case that seems to appear frequently for our investigation:
  # pt_reduced of the form a * x for non-zero number a and remaining variables x, y subject to a reduced SR generator x * y and scaling relation [*, != 0].
  # This only gives the solution [0:1], so one intersection point.
  if length(gs_reduced) == 1 && length(remaining_vars) == 2
    mons_list = collect(monomials(pt_reduced))
    if length(mons_list) == 1 && collect(coefficients(pt_reduced))[1] != 0
      list_of_exps = collect(exponents(mons_list[1]))[1]
      number_of_zeros = count(==(0), list_of_exps)
      if number_of_zeros == length(list_of_exps) - 1
        if gs_reduced[1] == remaining_vars[1] * remaining_vars[2]
          if (mons_list[1] == remaining_vars[1] && reduced_scaling_relations[1,2] != 0) || (mons_list[1] == remaining_vars[2] && reduced_scaling_relations[1,1] != 0)
            if list_of_exps[findfirst(x -> x > 0, list_of_exps)] == 1
              inter_dict[indices] = ZZ(1)
              return ZZ(1)
            end
          end
        end
      end
    end
  end
    
  # C.7 Cover a case that seems to appear frequently for our investigation. It looks as follows:
  # pt_reduced = -5700*w8*w10
  # remaining_vars = MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[w8, w10]
  # gs_reduced = MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[w8*w10]
  # reduced_scaling_relations = [1 1]
  # This gives exactly two solutions, namely [0:1] and [1:0].
  if length(gs_reduced) == 1 && length(remaining_vars) == 2
    mons_list = collect(monomials(pt_reduced))
    if length(mons_list) == 1 && mons_list[1] == remaining_vars[1] * remaining_vars[2]
      if gs_reduced[1] == remaining_vars[1] * remaining_vars[2]
        if reduced_scaling_relations[1,1] != 0 && reduced_scaling_relations[1,2] != 0
          inter_dict[indices] = 2
          return 2
        end
      end
    end
  end

  # C.8 Check if this was covered in our special cases
  if haskey(s_inter_dict, string([pt_reduced, gs_reduced, remaining_vars, reduced_scaling_relations]))
    numb = s_inter_dict[string([pt_reduced, gs_reduced, remaining_vars, reduced_scaling_relations])]
    inter_dict[indices] = numb
    return numb
  end

  # C.9 In all other cases, proceed via a rationally equivalent cycle
  #=
  println("")
  println("FOUND CASE THAT CANNOT YET BE DECIDED!")
  println("$pt_reduced")
  println("$remaining_vars")
  println("$gs_reduced")
  println("$indices")
  println("$reduced_scaling_relations")
  println("TRYING WITH EQUIVALENT CYCLE")
  println("")
  =#
  numb = intersection_from_equivalent_cycle(v, indices, hypersurface_equation, inter_dict, s_inter_dict, data)
  s_inter_dict[string([pt_reduced, gs_reduced, remaining_vars, reduced_scaling_relations])] = numb
  return numb

end



# ---------------------------------------------------------------------------------------------------------
# (2) Compute the intersection product from a rationally equivalent cycle.
# ---------------------------------------------------------------------------------------------------------

function intersection_from_equivalent_cycle(v::NormalToricVariety, indices::NTuple{4, Int64}, hypersurface_equation::MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}, inter_dict::Dict{NTuple{4, Int64}, ZZRingElem}, s_inter_dict::Dict{String, ZZRingElem}, data::NamedTuple)
  coeffs_list, tuple_list = Oscar._rationally_equivalent_cycle(v, indices, data)
  intersect_numb = 0
  for k in 1:length(tuple_list)
    if !is_zero(coeffs_list[k])
      intersect_numb += coeffs_list[k] * sophisticated_intersection_product(v, tuple_list[k], hypersurface_equation, inter_dict, s_inter_dict, data)
    end
  end
  @req is_integer(intersect_numb) "Should have expected to find only integer intersection numbers, but found $intersect_numb for $indices"
  inter_dict[indices] = ZZ(intersect_numb)
  return ZZ(intersect_numb)
end



# ---------------------------------------------------------------------------------------------------------
# (3) A function to reduce the hypersurface polynomial to {xi = 0} with i the entries of the tuple indices.
# ---------------------------------------------------------------------------------------------------------

function _reduce_hypersurface_equation(v::NormalToricVariety, p_hyper::MPolyRingElem, indices::NTuple{4, Int64}, data::NamedTuple)

  # Set variables to zero in the hypersurface equation
  vanishing_vars_pos = unique(indices)  
  new_p_hyper = divrem(p_hyper, data.gS[vanishing_vars_pos[1]])[2]
  for m in 2:length(vanishing_vars_pos)
    new_p_hyper = divrem(new_p_hyper, data.gS[vanishing_vars_pos[m]])[2]
  end

  # Is the resulting polynomial constant?
  if is_constant(new_p_hyper)
    return [new_p_hyper, [], [], zero_matrix(ZZ, 0, 0)]
  end

  # Identify the remaining variables
  remaining_vars_pos = Set(1:length(data.gS))
  for my_exps in data.sr_ideal_pos
    len_my_exps = length(my_exps)
    inter_len = count(idx -> idx in vanishing_vars_pos, my_exps)
    if len_my_exps == inter_len + 1
      delete!(remaining_vars_pos, my_exps[findfirst(idx -> !(idx in vanishing_vars_pos), my_exps)])
    end
  end
  set_to_one_list = sort([k for k in 1:length(data.gS) if k ∉ remaining_vars_pos])
  remaining_vars_pos = setdiff(collect(remaining_vars_pos), vanishing_vars_pos)
  remaining_vars = [data.gS[k] for k in remaining_vars_pos]

  # Extract remaining Stanley-Reisner ideal relations
  sr_reduced = Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}()
  for k in 1:length(data.sr_ideal_pos)
    if isdisjoint(set_to_one_list, data.sr_ideal_pos[k])
      new_gen = prod(data.gS[m] for m in data.sr_ideal_pos[k] if m ∉ vanishing_vars_pos)
      if !(new_gen in sr_reduced)
        push!(sr_reduced, new_gen)
      end
    end
  end

  # Identify the remaining scaling relations
  prepared_scaling_relations = zero_matrix(ZZ, length(data.scalings[1]), length(set_to_one_list) + length(remaining_vars_pos))
  for k in 1:(length(set_to_one_list) + length(remaining_vars_pos))
    col = k <= length(set_to_one_list) ? data.scalings[set_to_one_list[k]] : data.scalings[remaining_vars_pos[k - length(set_to_one_list)]]
    for l in 1:length(col)
      prepared_scaling_relations[l, k] = col[l]
    end
  end
  prepared_scaling_relations = hnf(prepared_scaling_relations)
  reduced_scaling_relations = prepared_scaling_relations[length(set_to_one_list) + 1: nrows(prepared_scaling_relations), length(set_to_one_list) + 1 : ncols(prepared_scaling_relations)]

  # Identify the final form of the reduced hypersurface equation, by setting all variables to one that we can
  images = [k in remaining_vars_pos ? data.gS[k] : one(data.S) for k in 1:length(data.gS)]
  pt_reduced = evaluate(new_p_hyper, images)

  # Return the result
  return [pt_reduced, sr_reduced, remaining_vars, reduced_scaling_relations]
end



# ---------------------------------------------------------------------------------------------------------
# (4) A function to find a rationally equivalent algebraic cycle.
# ---------------------------------------------------------------------------------------------------------

function _rationally_equivalent_cycle(v::NormalToricVariety, indices::NTuple{4, Int64}, data::NamedTuple)

  # Identify positions of the single and triple variable
  power_variable = indices[rand(1:length(indices))]
  for k in Set(indices)
    if count(==(k), indices) > 1
      power_variable = k
      break
    end
  end
  other_variables = collect(Set(filter(!=(power_variable), indices)))
  pos_power_variable = findfirst(==(power_variable), indices)

  # Let us simplify the problem by extracting the entries in the columns of single_variables and double_variables of the linear relation matrix
  simpler_matrix = matrix(QQ, data.linear_relations[[other_variables..., power_variable], :])
  b = zero_matrix(QQ, length(other_variables) + 1, 1)
  b[nrows(b), 1] = 1
  A = solve(simpler_matrix, b; side =:right)
  
  # Now form the relation in case...
  employed_relation = -sum((data.linear_relations[:, k] .* A[k,1]) for k in 1:nrows(A))
  employed_relation[power_variable] = 0

  # Populate `coeffs` and `tuples` that form rationally equivalent algebraic cycle
  coeffs = Vector{QQFieldElem}()
  tuples = Vector{NTuple{4, Int64}}()
  for k in 1:length(employed_relation)
    if employed_relation[k] != 0

      # Form the new tuple
      new_tuple = (indices[1:pos_power_variable-1]..., k, indices[pos_power_variable+1:end]...)
      new_tuple = Tuple(sort(collect(new_tuple)))

      pos = findfirst(==(new_tuple), tuples)
      if pos !== nothing
        coeffs[pos] += employed_relation[k]
      else
        push!(coeffs, employed_relation[k])
        push!(tuples, Tuple(new_tuple))
      end

    end
  end

  return [coeffs, tuples]

end
