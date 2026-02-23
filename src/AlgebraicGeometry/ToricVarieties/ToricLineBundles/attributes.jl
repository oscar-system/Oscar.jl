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
    sheaf_cohomology(l::ToricLineBundle, i::Int; algorithm::Symbol=:cohomcalg)

Compute the dimension of the *i*-th sheaf cohomology group of the toric line bundle `l`.

The keyword argument `algorithm::Symbol` selects the algorithm used for the computation.

Allowed values are:
  - `:cohomcalg` — use the cohomCalg algorithm (cf. [BJRR10](@cite), [BJRR10*1](@cite)).
    Requires the toric variety to be simplicial and projective.

  - `:chamber` — chamber counting algorithm (cf. [CLS11](@cite), p.398).
    Requires the toric variety to be simplicial and complete.

  - `:local` — local cohomology method (cf. [CLS11](@cite), Section 9.5).

By default, `:cohomcalg` is used.

# Examples
```jldoctest
julia> dP3 = del_pezzo_surface(NormalToricVariety, 3)
Normal toric variety

julia> sheaf_cohomology(toric_line_bundle(dP3, [4, 1, 1, 1]), 0)
12

julia> sheaf_cohomology(toric_line_bundle(dP3, [4, 1, 1, 1]), 0, algorithm = :chamber)
12

julia> sheaf_cohomology(toric_line_bundle(dP3, [4, 1, 1, 1]), 0, algorithm = :local)
12
```
"""
function sheaf_cohomology(l::ToricLineBundle, i::Int; algorithm::Symbol=:cohomcalg)
  v = toric_variety(l)
  has_attribute(v, :vanishing_sets) && contains(vanishing_sets(v)[i + 1], l) && return ZZ(0)
  return sheaf_cohomology(l; algorithm=algorithm)[i + 1]
end

@doc raw"""
    sheaf_cohomology(l::ToricLineBundle; algorithm::Symbol=:cohomcalg)

Compute the dimension of the sheaf cohomology groups of the toric line bundle `l`.

The keyword argument `algorithm::Symbol` selects the algorithm used for the computation.

Allowed values are:
  - `:cohomcalg` — use the cohomCalg algorithm (cf. [BJRR10](@cite), [BJRR10*1](@cite)).
    Requires the toric variety to be simplicial and projective.

  - `:chamber` — chamber counting algorithm (cf. [CLS11](@cite), p.398).
    Requires the toric variety to be simplicial and complete.

  - `:local` — local cohomology method (cf. [CLS11](@cite), Section 9.5).

By default, `:cohomcalg` is used.

# Examples
```jldoctest
julia> dP3 = del_pezzo_surface(NormalToricVariety, 3)
Normal toric variety

julia> sheaf_cohomology(toric_line_bundle(dP3, [1, 2, 3, 4]))
3-element Vector{ZZRingElem}:
 0
 16
 0

julia> sheaf_cohomology(toric_line_bundle(dP3, [1, 2, 3, 4]); algorithm = :chamber)
3-element Vector{ZZRingElem}:
 0
 16
 0

julia> sheaf_cohomology(toric_line_bundle(dP3, [1, 2, 3, 4]); algorithm = :local)
3-element Vector{ZZRingElem}:
 0
 16
 0

julia> sheaf_cohomology(toric_line_bundle(dP3, [-3,-2,-2,-2]))
3-element Vector{ZZRingElem}:
 0
 2
 0

julia> sheaf_cohomology(toric_line_bundle(dP3, [-3,-2,-2,-2]); algorithm = :chamber)
3-element Vector{ZZRingElem}:
 0
 2
 0

julia> sheaf_cohomology(toric_line_bundle(dP3, [-3,-2,-2,-2]); algorithm = :local)
3-element Vector{ZZRingElem}:
 0
 2
 0
```
"""
@attr Vector{ZZRingElem} function sheaf_cohomology(l::ToricLineBundle; algorithm::Symbol=:cohomcalg)
  v = toric_variety(l)
  if algorithm === :cohomcalg
    @req (is_simplicial(v) && is_projective(v)) "the currently implemented cohomCalg algorithm only applies to toric varieties that are simplicial and projective"
    return all_cohomologies_via_cohomcalg(l)
  elseif algorithm === :chamber
    @req (is_complete(v) && is_simplicial(v)) "the chamber counting algorithm only applies to toric varieties that are simplicial and complete"
    return _all_cohomologies_via_cech(l)
  elseif algorithm === :local
    ctx = local_cohomology_context_object(v)
    d = divisor_class(toric_divisor_class(l))
    coh = cohomology_model(ctx, d)
    return ZZRingElem[ZZ(ngens(coh[i])) for i in 0:-1:(-dim(v))]
  else
    throw(ArgumentError("Unknown algorithm=$algorithm. Use :cohomcalg, :chamber, or :local."))
  end
end
