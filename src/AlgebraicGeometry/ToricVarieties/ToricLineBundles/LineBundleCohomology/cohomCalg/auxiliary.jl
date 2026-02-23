###################################
# (0) Prepare to cohomCalg to compute line bundle cohomolgy
###################################

function all_cohomologies_via_cohomcalg(l::ToricLineBundle)::Vector{ZZRingElem}
  v = toric_variety(l)
  @req (is_simplicial(v) && is_projective(v)) "the currently implemented cohomCalg algorithm only applies to toric varieties that are simplicial and projective"
  class = vec([ZZRingElem(x) for x in divisor_class(toric_divisor_class(l)).coeff])
  
  # Minimal example:
  #
  # The following lines provide an interface to cohomCalg. This is easiest explained with an example.
  # 
  # Example: Compute the cohomologies of the structure sheaf of the 2-dimensional projective space.
  #
  # Step 1: We trigger the following command so that cohomCalg computes these cohomologies for us:
  #
  # ./cohomCalg "vertex x1 | GLSM: ( 1 ); vertex x2 | GLSM: ( 1 ); vertex x3 | GLSM: ( 1 ); srideal [x1*x2*x3]; ambientcohom O(0 );"
  #
  # The notation GLSM is for generalized linear sigma model.
  # This is a physics model which famously triggered interest in toric geometry for physics applications.
  # The cohomCalg algorithm relies only on the grading of the Cox ring and the Stanley-Reisner ideal.
  # <-> Those two pieces of information are (more or less) the defining data of a GLSM.
  #
  # Step 2: We receive the following return value from cohomCalg:
  #
  # "{True, {{1,0,0}, {{0, 1*1}}}}"
  #
  # The "True" tells us that the run was successful. Otherwise, we would find "False".
  # 
  # Whenever we encounter "False", this most likely indicates wrong input data. This could be any of the following:
  # (a) wrong variable names: cohomCalg accepts "x1, x2, ..." (and also "y1, y2, ..."), but "_x1, _x2" or "x[1]", "x[2]" will cause errors.
  # (b) an inconsistent grading of the Cox ring: This could e.g. lead to infinitely many monoms of a fixed degree.
  # (c) an inconsistent Stanley-Reisner ideal.
  # 
  # For sake of simplicity, our implementation creates a dictionary that maps "our" variable names to "x1, x2, ...".
  # The latter are then passed to cohomCalg.
  # 
  # After this boolean, the above return value contains the line bundle cohomologies and some intermediate results.
  # The first argument in the { }-brackets are the line bundle cohomologies in question. Here: {1,0,0}.
  # 
  # By the cohomCalg algorithm, the line bundle cohomologies are given by certain rationoms, that is fractions of 
  # monoms in the Cox ring of the variety. Behind the scenes, cohomCalg reads-off monomials which contribute
  # to the denominators of these rationoms. For sake of simplicity, let us refer to these monomials as "partial-denominators".
  # 
  # Crucially: These partial-denominators are read-off from the Stanley-Reisner ideal. Hence, they encode information about
  # the toric space and are not specific to a certain line bundle. In https://arxiv.org/abs/1802.08860, this is used to infer
  # refined vanishing sets of line bundle cohomologies on toric spaces. Implementations of which are currently available in
  # https://github.com/homalg-project/ToricVarieties_project and will be migrated to OSCAR soon.
  # 
  # To make contact with a line bundle in question, pick a partial-denominator and look for all rationoms with multi-degree
  # matching that of the line bundle in question. In particular, under the assumptions that the variety in question be either
  # smooth, complete or alternatively simplicial, projective, one finds a finite number of such rationoms for each partial-denominator.
  # 
  # Once these rationoms are identified, two questions remain:
  # (1) To which line bundle cohomology do they contribute?
  # (2) Could they encode more than a one-dimensional linear subspace of these cohomology spaces?
  # 
  # The answers to both questions are involved and described in detail in the original work https://arxiv.org/pdf/1003.5217.pdf
  # and the proofs of this algorithm in https://arxiv.org/abs/1006.2392, https://arxiv.org/abs/1006.0780. In particular note that
  # the answer to (2) involved the evaluation of certain "remnant cohomologies", so that a single rationom may indeed encode
  # a vector space of line bundle cohomologies which is strictly larger than 1.
  # 
  # Let us now return to the above return value {{0, 1*1}}. The fact that this list contains only a single sub-list states that only
  # a single partial-denominator was identified. The first value tells to which cohomology group it contributes. Here we find 0,
  # so it contributes to H^0. The second value is, in the following order, the product of the multiplicity (as in response to 
  # (2) above) and the exact form of the partial-denominator. So here the multiplicity is 1 and the partial-denominator is 1.
  # For other examples, we could for example find "2 * x2 * x3 * x5", which would mean that the multiplicity is 2
  # and the partial denominator in question is "x2 * x3 * x5".
  # 
  # For more details, please refer to the cohomCalg manual: https://arxiv.org/pdf/1003.5217.pdf.
  # 
  # Step 3: Extract "True" (or "False") and the first argument of the following "{ }" by suitable parsing.
  #
  # -> Hooray! We found the line bundle cohomologies in question.

  # obtain the command string
  class = vec([ZZRingElem(x) for x in divisor_class(toric_divisor_class(l)).coeff])
  command = command_string(v, class)

  # execute cohomCalg
  out = Pipe()
  err = Pipe()
  process = run(
    pipeline(
      ignorestatus(`$(cohomCalg_jll.cohomcalg()) --integrated --in=$(command)`);
      stdout=out,
      stderr=err,
    ),
  )
  close(out.in)
  close(err.in)

  # was there an error?
  stderr = read(err, String)
  code = process.exitcode
  if code != 0
    error("cohomCalg encountered the error " * stderr)
  end

  # read out the result
  stdout = read(out, String)
  result = [parse(ZZRingElem, c) for c in split(chop(chop(split(stdout, "{")[4])), ",")]

  # consistency check
  if length(result) != dim(v) + 1
    error(
      "cohomCalg should return list of length $(dim(v)+1) but returned list of length $(length(result))"
    )
  end

  # return result
  return result

end

###################################
# (1) Call cohomCalg to compute line bundle cohomolgy
###################################

function command_string(v::NormalToricVarietyType, c::Vector{ZZRingElem})
  # Initialize a list of strings, which we will eventually join to form the command string
  string_list = Vector{String}()

  # Define helper function
  joincomma(list) = join([string(x) for x in list], ", ")

  # Add information about grading of Cox ring to string_list
  divisors = gens(torusinvariant_weil_divisor_group(v))
  for i in 1:length(divisors)
    tmp = joincomma(
      map_from_torusinvariant_weil_divisor_group_to_class_group(v)(divisors[i]).coeff
    )
    push!(string_list, "vertex x$i|GLSM:($(tmp))")
  end

  # Add information about the Stanley-Reisner ideal to string_list
  current_coordinate_names = [string(x) for x in Hecke.gens(cox_ring(v))]
  new_coordinate_names = [Symbol(:x, i) for i in 1:length(current_coordinate_names)]
  new_ring, _ = polynomial_ring(coefficient_ring(v), new_coordinate_names; cached=false)
  generators = [string(g) for g in gens(stanley_reisner_ideal(new_ring, v))]
  push!(string_list, "srideal [" * joincomma(generators) * "]")

  # Add line bundle information to string_list
  push!(string_list, "ambientcohom O(" * joincomma(c) * ");")

  # Join and return
  return join(string_list, ";")
end
command_string(v::NormalToricVarietyType) = command_string(
  v, [ZZRingElem(0) for i in 1:torsion_free_rank(picard_group_with_map(v)[1])]
)

###################################
# (2) Compute the denominators of the rationoms
###################################

function contributing_denominators(variety::NormalToricVarietyType)

  # execute cohomCalg
  out = Pipe()
  err = Pipe()
  process = run(
    pipeline(
      ignorestatus(
        `$(cohomCalg_jll.cohomcalg()) --verbose2 --in=$(command_string(variety))`
      );
      stdout=out,
      stderr=err,
    ),
  )
  close(out.in)
  close(err.in)

  # error during execution?
  stderr = read(err, String)
  code = process.exitcode
  if code != 0
    error("cohomCalg encountered the error " * stderr)
  end

  # process result
  stdout = read(out, String)
  start = findfirst("Final list of contributing monomials with factors:", stdout)[1]
  diff = findfirst("Verbose Level 1:", stdout)[1]
  output_string_reduced = [strip(s) for s in split(SubString(stdout, start, diff), "\n")]
  output_string_reduced = output_string_reduced[4:(length(output_string_reduced) - 3)]

  # ambiguous monomial contributions found during execution?
  if last(output_string_reduced) !=
    "There are no ambiguous contribution monomials to the variety."
    throw(ArgumentError("cohomCalg encountered ambiguous monomial contributions"))
  end

  # Monomials can in principal come with contributing factor different from 1.
  # In this case, the word "factor ..." occurs. Here, we do not care about these factors.
  # To process the result efficiently, we remove this information.
  for i in 1:length(output_string_reduced)
    if occursin("factor", output_string_reduced[i])
      output_string_reduced[i] = strip(
        SubString(
          output_string_reduced[i], 1, findfirst("factor", output_string_reduced[i])[1] - 1
        ),
      )
    end
  end

  # at what positions in the output_string does the list of contributions to cohomology group i end and that of i+1 begins?
  hlist = Int[]
  for i in 1:(length(output_string_reduced) - 1)
    if output_string_reduced[i][1] == 'H'
      push!(hlist, i)
    end
  end

  # form list of the contributing monomials
  contributing_monomials = [[""] for i in 1:(dim(variety) + 1)]
  for i in 1:(dim(variety) + 1)
    if ((i <= dim(variety)) && (hlist[i] + 1 !== hlist[i + 1]))
      contributing_monomials[i] = [
        output_string_reduced[j] for j in (hlist[i] + 1):(hlist[i + 1] - 1)
      ]
    elseif ((i == dim(variety) + 1) && (hlist[i] + 1 !== length(output_string_reduced)))
      contributing_monomials[i] = [
        output_string_reduced[j] for j in (hlist[i] + 1):(length(output_string_reduced) - 1)
      ]
    end
  end

  # we execute cohomCalg with homogeneous variable names "x1", "x2" and so on
  # -> translate back into actually chosen variable names
  gens = Hecke.gens(cox_ring(variety))
  for cm in contributing_monomials
    if cm == [""]
      empty!(cm)
    else
      for j in 1:length(cm)
        m = cm[j]
        present_variables = [occursin("x$k", m) for k in 1:ngens(cox_ring(variety))]
        cm[j] = string(
          prod(gens[k]^present_variables[k] for k in 1:ngens(cox_ring(variety)))
        )
      end
    end
  end

  # return result
  return contributing_monomials
end

######################################
# (3) Turn denominators into shifted cones (-> vanishing sets)
######################################

function turn_denominator_into_polyhedron(variety::NormalToricVarietyType, monom::String)

  # (1) which variables appear in the monom?
  present_variables = [
    occursin(string(coordinate_names(variety)[i]), monom) for
    i in 1:ngens(cox_ring(variety))
  ]

  # (2) compute generators of the semigroup
  weights = [k.coeff for k in Oscar._cox_ring_weights(variety)]
  gens = reduce(
    vcat,
    unique([
      (-1)^Int(present_variables[i]) * weights[i] for i in 1:length(present_variables)
    ]),
  )

  # (3) compute offset
  offset = zero(parent(weights[1]))
  for i in 1:length(present_variables)
    if present_variables[i]
      offset -= weights[i]
    end
  end
  return convex_hull(offset, gens)
end
