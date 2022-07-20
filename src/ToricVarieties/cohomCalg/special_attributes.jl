###########################
# (1) Special attributes of toric varieties
###########################

@doc Markdown.doc"""
    vanishing_sets(variety::AbstractNormalToricVariety)

Compute the vanishing sets of an abstract toric variety `v` by use of the cohomCalg algorithm.
"""
@attr function vanishing_sets(variety::AbstractNormalToricVariety)
    denominator_contributions = contributing_denominators(variety)
    vs = ToricVanishingSet[]
    for i in 1:length(denominator_contributions)
        list_of_polyhedra = [turn_denominator_into_polyhedron(variety, m) for m in denominator_contributions[i]]
        push!(vs, ToricVanishingSet(variety, list_of_polyhedra, i-1))
    end
    return vs::Vector{ToricVanishingSet}
end
export vanishing_sets


###########################
# (2) Special attributes of toric line bundles
###########################

@doc Markdown.doc"""
    all_cohomologies(l::ToricLineBundle)

Computes the dimension of all sheaf cohomologies of the 
toric line bundle `l` by use of the cohomCalg algorithm 
[BJRR10](@cite),
[cohomCalg:Implementation(@cite),
[RR10](@cite),
[Jow11](@cite),
[BJRR12](@cite).

# Examples
```jldoctest
julia> dP3 = del_pezzo(3)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> all_cohomologies(ToricLineBundle(dP3, [1,2,3,4]))
3-element Vector{fmpz}:
 0
 16
 0
```
"""
function all_cohomologies(l::ToricLineBundle)
    # check if we can apply cohomCalg
    v = toric_variety(l)
    if !((is_smooth(v) && is_complete(v)) || (is_simplicial(v) && is_projective(v)))
        throw(ArgumentError("cohomCalg only applies to toric varieties that are either smooth, complete or simplicial, projective"))
    end
    
    # extract vector of currently-known cohomology dimensions (or create it if necessary)
    return get_attribute!(l, :all_cohomologies) do
        
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
        # "{True,{{1,0,0},{{0,1*1}}}}"
        #
        # The "True" tells us that the run was successful. Otherwise, we would find "False".
        # 
        # Whenever we encounter "False", this most likely indicates wrong input data. This could be any of the following:
        # (a) wrong variable names: cohomCalg accepts "x1, x2, ..." (and also "y1, y2, ..."), but "_x1, _x2" or "x[1]", "x[2]" will cause errors.
        # (b) an inconsistent grading of the Cox ring: This could e.g. lead to infinitely many monoms of a fixed degree.
        # (c) an inconsistent Stanley-Reisner ideal.
        # 
        # For sake of simplicity, our implementation creates a dictionary that maps "our" variable names to "x1, x2,...".
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
        # Let us now return to the above return value {{0,1*1}}. The fact that this list contains only a single sub-list states that only 
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
        class = vec([fmpz(x) for x in divisor_class(l).coeff])
        command = command_string(v, class)
        
        # execute cohomCalg
        out = Pipe()
        err = Pipe()
        process = run(pipeline(ignorestatus(`$(cohomCalg_jll.cohomcalg()) --integrated --in=$(command)`), stdout=out, stderr=err))
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
        result = [fmpz(parse(Int,c)) for c in split(chop(chop(split(stdout, "{" )[4])), ",")]
        
        # consistency check
        if length(result) != dim(v)+1
            error("cohomCalg should return list of length $(dim(v)+1) but returned list of length $(length(result))")
        end
        
        # return result
        return result
    end
end
export all_cohomologies


@doc Markdown.doc"""
    cohomology(l::ToricLineBundle, i::Int)

Computes the dimension of the i-th sheaf cohomology of the
toric line bundle `l` by use of the cohomCalg algorithm
[BJRR10](@cite),
[cohomCalg:Implementation(@cite),
[RR10](@cite),
[Jow11](@cite),
[BJRR12](@cite).

# Examples
```jldoctest
julia> dP3 = del_pezzo(3)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> cohomology(ToricLineBundle(dP3, [4,1,1,1]), 0)
12
```
"""
function cohomology(l::ToricLineBundle, i::Int)
    v = toric_variety(l)
    if has_attribute(v, :vanishing_sets)
        tvs = vanishing_sets(v)[i+1]
        if contains(tvs, l)
            return 0
        end
    end
    return all_cohomologies(l)[i+1]
end
export cohomology
