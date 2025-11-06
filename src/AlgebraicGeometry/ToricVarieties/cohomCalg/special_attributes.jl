###########################
# (1) Special attributes of toric varieties
###########################R

@doc raw"""
    vanishing_sets(variety::NormalToricVarietyType)

Compute the vanishing sets of an abstract toric variety `v` by use of the cohomCalg algorithm.
"""
@attr Vector{ToricVanishingSet} function vanishing_sets(variety::NormalToricVarietyType)
  @req (is_simplicial(variety) && is_projective(variety)) "the currently implemented cohomCalg algorithm only applies to toric varieties that are simplicial and projective"
  denominator_contributions = contributing_denominators(variety)
  vs = ToricVanishingSet[]
  for i in 1:length(denominator_contributions)
    list_of_polyhedra = Polyhedron{QQFieldElem}[
      turn_denominator_into_polyhedron(variety, m) for m in denominator_contributions[i]
    ]
    push!(vs, ToricVanishingSet(variety, list_of_polyhedra, [i - 1]))
  end
  return vs
end

@doc raw"""
    immaculate_line_bundles(variety::NormalToricVarietyType)

Compute all immaculate line bundles as a toric vanishing set by
intersecting the vanishing sets for all cohomology indices.

# Examples
```jldoctest
julia> dP1 = del_pezzo_surface(NormalToricVariety, 1)
Normal toric variety

julia> ilb = immaculate_line_bundles(dP1)
Toric vanishing set for cohomology indices [0, 1, 2]

julia> polyhedra(ilb)
4-element Vector{Polyhedron{QQFieldElem}}:
 Polyhedron in ambient dimension 2
 Polyhedron in ambient dimension 2
 Polyhedron in ambient dimension 2
 Polyhedron in ambient dimension 2

julia> print_constraints(polyhedra(ilb)[1])
-x_1 <= 0
-x_1 + x_2 <= 0

julia> print_constraints(polyhedra(ilb)[2])
-x_1 + x_2 <= 0
x_2 <= -2

julia> print_constraints(polyhedra(ilb)[3])
-x_2 <= -1
x_1 - x_2 <= -2

julia> print_constraints(polyhedra(ilb)[4])
x_1 - x_2 <= -2
x_1 <= -3
```
"""
@attr Any function immaculate_line_bundles(variety::NormalToricVarietyType)
  @req (is_simplicial(variety) && is_projective(variety)) "the currently implemented cohomCalg algorithm only applies to toric varieties that are simplicial and projective"
  denominator_contributions = reduce(vcat, contributing_denominators(variety))
  list_of_polyhedra = Polyhedron{QQFieldElem}[
    turn_denominator_into_polyhedron(variety, m) for m in denominator_contributions
  ]
  return ToricVanishingSet(variety, list_of_polyhedra, collect(0:dim(variety)))
end

###########################
# (2) Special attributes of toric line bundles
###########################

@doc raw"""
    sheaf_cohomology(l::ToricLineBundle)

Compute the dimension of all sheaf cohomologies of the 
toric line bundle `l`. The default algorithm is the cohomCalg algorithm 
[BJRR10](@cite), [BJRR10*1](@cite) (see also [RR10](@cite),
[Jow11](@cite) and [BJRR12](@cite)). It is also possible to specify algorithm = "chamber counting"
in which case the chamber counting algorithm will be used [CLS11](@cite) p.398, 
or to use local cohomology as in [CLS11](@cite), Section 9.5.

# Examples
```jldoctest
julia> dP3 = del_pezzo_surface(NormalToricVariety, 3)
Normal toric variety

julia> sheaf_cohomology(toric_line_bundle(dP3, [1, 2, 3, 4]))
3-element Vector{ZZRingElem}:
 0
 16
 0

julia> sheaf_cohomology(toric_line_bundle(dP3, [1, 2, 3, 4]); algorithm = "chamber counting")
3-element Vector{ZZRingElem}:
 0
 16
 0

julia> sheaf_cohomology(toric_line_bundle(dP3, [-3,-2,-2,-2]))
3-element Vector{ZZRingElem}:
 0
 2
 0

julia> sheaf_cohomology(toric_line_bundle(dP3, [-3,-2,-2,-2]); algorithm = "chamber counting")
3-element Vector{ZZRingElem}:
 0
 2
 0
```
"""
@attr Vector{ZZRingElem} function sheaf_cohomology(
  l::ToricLineBundle; algorithm::String="cohomCalg"
)
  v = toric_variety(l)
  if occursin("cohomcalg", lowercase(algorithm))
    # check if we can apply cohomCalg
    @req (is_simplicial(v) && is_projective(v)) "the currently implemented cohomCalg algorithm only applies to toric varieties that are simplicial and projective"

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
  elseif occursin("chamber", lowercase(algorithm))
    @req (is_complete(v) && is_simplicial(v)) "the chamber counting algorithm only applies to toric varieties that are simplicial and complete"
    return _all_cohomologies_via_cech(l)
  elseif occursin("local", lowercase(algorithm))
    ctx = local_cohomology_context_object(v)
    d = divisor_class(toric_divisor_class(l))
    coh = cohomology_model(ctx, d)
    return ZZRingElem[ZZ(ngens(coh[i])) for i in 0:-1:-dim(v)]
  end
end

@doc raw"""
    sheaf_cohomology(l::ToricLineBundle, i::Int)

Compute the dimension of the i-th sheaf cohomology of the
toric line bundle `l` by use of the cohomCalg algorithm
[BJRR10](@cite), [BJRR10*1](@cite) (see also [RR10](@cite),
[Jow11](@cite) and [BJRR12](@cite)).

# Examples
```jldoctest
julia> dP3 = del_pezzo_surface(NormalToricVariety, 3)
Normal toric variety

julia> sheaf_cohomology(toric_line_bundle(dP3, [4, 1, 1, 1]), 0)
12
```
"""
function sheaf_cohomology(l::ToricLineBundle, i::Int; algorithm::String="cohomCalg")
  has_attribute(l, :sheaf_cohomology) && return get_attribute(l, :sheaf_cohomology)::Vector{ZZRingElem}[i]
  table = get_attribute!(l, :sheaf_cohomology_table) do
    Dict{Int, ZZRingElem}()
  end
  return get!(table, i) do
    _sheaf_cohomology(l, i; algorithm)
  end
end

function _sheaf_cohomology(l::ToricLineBundle, i::Int; algorithm::String="cohomCalg")
  v = toric_variety(l)
  if has_attribute(v, :vanishing_sets)
    tvs = vanishing_sets(v)[i + 1]
    if contains(tvs, l)
      return 0
    end
  end
  if occursin("cohomcalg", lowercase(algorithm))
    return sheaf_cohomology(l; algorithm="cohomCalg")[i + 1]
  elseif occursin("chamber", lowercase(algorithm))
    return sheaf_cohomology(l; algorithm="chamber")[i + 1]
  elseif occursin("local", lowercase(algorithm))
    ctx = local_cohomology_context_object(v)
    d = divisor_class(toric_divisor_class(l))
    coh = cohomology_model(ctx, d)
    return ZZ(ngens(coh[-i]))
  end
end

