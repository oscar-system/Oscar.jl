#####################
# 1. Defining data of a line bundle
#####################

@doc raw"""
    picard_class(l::ToricLineBundle)

Return the class in the Picard group which defines the toric line bundle `l`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l = toric_line_bundle(v, [ZZRingElem(2)])
Toric line bundle on a normal toric variety

julia> picard_class(l)
Abelian group element [2]
```
"""
picard_class(l::ToricLineBundle) = l.picard_class

@doc raw"""
    toric_variety(l::ToricLineBundle)

Return the toric variety over which the toric line bundle `l` is defined.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l = toric_line_bundle(v, [ZZRingElem(2)])
Toric line bundle on a normal toric variety

julia> toric_variety(l)
Normal toric variety without torusfactor
```
"""
toric_variety(l::ToricLineBundle) = l.toric_variety

@doc raw"""
    toric_divisor(l::ToricLineBundle)

Return a toric divisor corresponding to the toric line bundle `l`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l = toric_line_bundle(v, [ZZRingElem(2)])
Toric line bundle on a normal toric variety

julia> toric_divisor(l)
Torus-invariant, cartier, non-prime divisor on a normal toric variety

julia> is_cartier(toric_divisor(l))
true
```
"""
@attr ToricDivisor function toric_divisor(l::ToricLineBundle)
  class = picard_class(l)
  map1 = map_from_torusinvariant_cartier_divisor_group_to_picard_group(toric_variety(l))
  map2 = map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(
    toric_variety(l)
  )
  image = map2(preimage(map1, class)).coeff
  coeffs = vec([ZZRingElem(x) for x in image])
  td = toric_divisor(toric_variety(l), coeffs)
  set_attribute!(td, :is_cartier, true)
  return td
end

@doc raw"""
    toric_divisor_class(l::ToricLineBundle)

Return a divisor class in the Class group corresponding to the toric line bundle `l`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l = toric_line_bundle(v, [ZZRingElem(2)])
Toric line bundle on a normal toric variety

julia> toric_divisor(l)
Torus-invariant, cartier, non-prime divisor on a normal toric variety

julia> is_cartier(toric_divisor(l))
true
```
"""
@attr ToricDivisorClass toric_divisor_class(l::ToricLineBundle) = toric_divisor_class(
  toric_divisor(l)
)

@doc raw"""
    degree(l::ToricLineBundle)

Return the degree of the toric line bundle `l`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l = toric_line_bundle(v, [ZZRingElem(2)])
Toric line bundle on a normal toric variety

julia> degree(l)
2
```
"""
@attr ZZRingElem degree(l::ToricLineBundle) = sum(coefficients(toric_divisor(l)))

#############################
# 2. Basis of global sections
#############################

@doc raw"""
    basis_of_global_sections_via_rational_functions(l::ToricLineBundle)

Return a basis of the global sections of the toric line bundle `l` in terms of rational functions.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l = toric_line_bundle(v, [ZZRingElem(2)])
Toric line bundle on a normal toric variety

julia> basis_of_global_sections_via_rational_functions(l)
6-element Vector{MPolyQuoRingElem{QQMPolyRingElem}}:
 x1_^2
 x2*x1_^2
 x2^2*x1_^2
 x1_
 x2*x1_
 1
```
"""
@attr Vector{MPolyQuoRingElem{QQMPolyRingElem}} function basis_of_global_sections_via_rational_functions(
  l::ToricLineBundle
)
  if has_attribute(toric_variety(l), :vanishing_sets)
    tvs = vanishing_sets(toric_variety(l))[1]
    if contains(tvs, l)
      return MPolyQuoRingElem{QQMPolyRingElem}[]
    end
  end
  characters = matrix(ZZ, lattice_points(polyhedron(toric_divisor(l))))
  return MPolyQuoRingElem{QQMPolyRingElem}[
    character_to_rational_function(
      toric_variety(l), vec([ZZRingElem(c) for c in characters[i, :]])
    ) for i in 1:nrows(characters)
  ]
end

@doc raw"""
    basis_of_global_sections_via_homogeneous_component(l::ToricLineBundle)

Return a basis of the global sections of the toric line bundle `l`
in terms of a homogeneous component of the Cox ring of `toric_variety(l)`.
For convenience, this method can also be called via
`basis_of_global_sections(l::ToricLineBundle)`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l = toric_line_bundle(v, [ZZRingElem(2)])
Toric line bundle on a normal toric variety

julia> basis_of_global_sections_via_homogeneous_component(l)
6-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x3^2
 x2*x3
 x2^2
 x1*x3
 x1*x2
 x1^2

julia> basis_of_global_sections(l)
6-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x3^2
 x2*x3
 x2^2
 x1*x3
 x1*x2
 x1^2
```
"""
@attr Vector{MPolyDecRingElem{QQFieldElem,QQMPolyRingElem}} function basis_of_global_sections_via_homogeneous_component(
  l::ToricLineBundle
)
  if has_attribute(toric_variety(l), :vanishing_sets)
    tvs = vanishing_sets(toric_variety(l))[1]
    if contains(tvs, l)
      return MPolyDecRingElem{QQFieldElem,QQMPolyRingElem}[]
    end
  end
  return monomial_basis(cox_ring(toric_variety(l)), divisor_class(toric_divisor_class(l)))
end
basis_of_global_sections(l::ToricLineBundle) =
  basis_of_global_sections_via_homogeneous_component(
    l
  )

#############################
# 3. Generic section
#############################

@doc raw"""
    generic_section(l::ToricLineBundle; range::UnitRange{Int64} = -10000:10000, rng::AbstractRNG = Random.default_rng())

Return a generic section of the toric line bundle `l`, that
is return the sum of all elements `basis_of_global_sections(l)`,
each multiplied by a random integer.

The optional keyword argument `range` can be used to set the range
of the random integers, e.g., `generic_section(l, range = -100:100)`

The random source used to create random coefficients can be set with
the optional argument `rng`.

# Examples
```jldoctest
julia> using Random;

julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l = toric_line_bundle(v, [ZZRingElem(2)])
Toric line bundle on a normal toric variety

julia> s = generic_section(l, rng = Random.Xoshiro(1234));

julia> parent(s) == cox_ring(toric_variety(l))
true
```
"""
function generic_section(
  l::ToricLineBundle;
  range::UnitRange{Int64}=-10000:10000,
  rng::AbstractRNG=Random.default_rng(),
)
  global_sections = basis_of_global_sections(l)

  if length(global_sections) == 0
    return zero(cox_ring(toric_variety(l)))
  end

  return sum(rand(rng, range) * b for b in global_sections)
end

###########################
# (4) Sheaf Cohomology
###########################

@doc raw"""
    sheaf_cohomology(l::ToricLineBundle, i::Int; algorithm::String="cohomCalg")

Compute the dimension of the i-th sheaf cohomology of the
toric line bundle `l`. The third argument allows to choose an algorithm.
By default, we employ the cohomCalg algorithm [BJRR10](@cite), [BJRR10*1](@cite),
see also [RR10](@cite), [Jow11](@cite) and [BJRR12](@cite).

It is also possible to specify algorithm = "chamber" in which case the chamber counting
algorithm will be used [CLS11](@cite) p.398. Also, it is possible to use local cohomology
as in [CLS11](@cite), Section 9.5 by specifying algorithm = "local".

# Examples
```jldoctest
julia> dP3 = del_pezzo_surface(NormalToricVariety, 3)
Normal toric variety

julia> sheaf_cohomology(toric_line_bundle(dP3, [4, 1, 1, 1]), 0, algorithm = "cohomCalg")
12

julia> sheaf_cohomology(toric_line_bundle(dP3, [4, 1, 1, 1]), 0, algorithm = "chamber")
12

julia> sheaf_cohomology(toric_line_bundle(dP3, [4, 1, 1, 1]), 0, algorithm = "local")
12
```
"""
function sheaf_cohomology(l::ToricLineBundle, i::Int; algorithm::String="cohomCalg")
  if has_attribute(l, :sheaf_cohomology)
    return ((get_attribute(l, :sheaf_cohomology)::Vector{ZZRingElem})[i + 1])::ZZRingElem
  end
  table = get_attribute!(l, :sheaf_cohomology_table) do
    Dict{Int,ZZRingElem}()
  end
  return get!(table, i) do
    v = toric_variety(l)
    if has_attribute(v, :vanishing_sets)
      contains(vanishing_sets(v)[i + 1], l) && return ZZ(0)
    end
    return sheaf_cohomology(l; algorithm)[i + 1]
  end
end

@doc raw"""
    sheaf_cohomology(l::ToricLineBundle; algorithm::String="cohomCalg")

Compute the dimension of all sheaf cohomologies of the 
toric line bundle `l`. The default algorithm is the cohomCalg algorithm 
[BJRR10](@cite), [BJRR10*1](@cite) (see also [RR10](@cite),
[Jow11](@cite) and [BJRR12](@cite)). 

It is also possible to specify algorithm = "chamber" in which case the chamber counting
algorithm will be used [CLS11](@cite) p.398. Also, it is possible to use local cohomology
as in [CLS11](@cite), Section 9.5 by specifying algorithm = "local".

# Examples
```jldoctest
julia> dP3 = del_pezzo_surface(NormalToricVariety, 3)
Normal toric variety

julia> sheaf_cohomology(toric_line_bundle(dP3, [1, 2, 3, 4]))
3-element Vector{ZZRingElem}:
 0
 16
 0

julia> sheaf_cohomology(toric_line_bundle(dP3, [1, 2, 3, 4]); algorithm = "chamber")
3-element Vector{ZZRingElem}:
 0
 16
 0

julia> sheaf_cohomology(toric_line_bundle(dP3, [1, 2, 3, 4]); algorithm = "local")
3-element Vector{ZZRingElem}:
 0
 16
 0

julia> sheaf_cohomology(toric_line_bundle(dP3, [-3,-2,-2,-2]))
3-element Vector{ZZRingElem}:
 0
 2
 0

julia> sheaf_cohomology(toric_line_bundle(dP3, [-3,-2,-2,-2]); algorithm = "chamber")
3-element Vector{ZZRingElem}:
 0
 2
 0

julia> sheaf_cohomology(toric_line_bundle(dP3, [-3,-2,-2,-2]); algorithm = "local")
3-element Vector{ZZRingElem}:
 0
 2
 0
```
"""
@attr Vector{ZZRingElem} function sheaf_cohomology(l::ToricLineBundle; algorithm::String="cohomCalg")
  v = toric_variety(l)
  if has_attribute(v, :sheaf_cohomology_table)
    table = get_attribute(v, :sheaf_cohomology_table)::Dict{Int,ZZRingElem}
    all(haskey(table, i) for i in 0:dim(v)) && return ZZRingElem[table[i] for i in 0:dim(v)]
  end
  if occursin("cohomcalg", lowercase(algorithm))
    @req (is_simplicial(v) && is_projective(v)) "the currently implemented cohomCalg algorithm only applies to toric varieties that are simplicial and projective"
    return _all_cohomologies_via_cohomcalg(l)
  elseif occursin("chamber", lowercase(algorithm))
    @req (is_complete(v) && is_simplicial(v)) "the chamber counting algorithm only applies to toric varieties that are simplicial and complete"
    return _all_cohomologies_via_cech(l)
  elseif occursin("local", lowercase(algorithm))
    ctx = local_cohomology_context_object(v)
    d = divisor_class(toric_divisor_class(l))
    coh = cohomology_model(ctx, d)
    return ZZRingElem[ZZ(ngens(coh[i])) for i in 0:-1:(-dim(v))]
  end
end

function _all_cohomologies_via_cohomcalg(l::ToricLineBundle)::Vector{ZZRingElem}
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
