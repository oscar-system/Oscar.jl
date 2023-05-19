################################################
# 1: Constructors with toric variety as base
################################################

@doc raw"""
    weierstrass_model(base::AbstractNormalToricVariety; completeness_check::Bool = true)

This method constructs a Weierstrass model over a given toric base
3-fold. The Weierstrass sections ``f`` and ``g`` are taken with (pseudo)random
coefficients.

# Examples
```jldoctest
julia> w = weierstrass_model(sample_toric_variety(); completeness_check = false)
Weierstrass model over a concrete base
```
"""
function weierstrass_model(base::AbstractNormalToricVariety; completeness_check::Bool = true)
  (f, g) = _weierstrass_sections(base)
  return weierstrass_model(f, g, base; completeness_check = completeness_check)
end


@doc raw"""
    weierstrass_model(f::MPolyRingElem, g::MPolyRingElem, base::AbstractNormalToricVariety; completeness_check::Bool = true)

This method operates analogously to `weierstrass_model(base::AbstractNormalToricVariety)`.
The only difference is that the Weierstrass sections ``f`` and ``g`` can be specified with non-generic values.

# Examples
```jldoctest
julia> base = sample_toric_variety()
Normal toric variety

julia> f = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^4)]);

julia> g = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)]);

julia> w = weierstrass_model(f, g, base; completeness_check = false)
Weierstrass model over a concrete base
```
"""
function weierstrass_model(f::MPolyRingElem, g::MPolyRingElem, base::AbstractNormalToricVariety; completeness_check::Bool = true)
  @req ((parent(f) == cox_ring(base)) && (parent(g) == cox_ring(base))) "All Weierstrass sections must reside in the Cox ring of the base toric variety"
  
  gens_base_names = [string(g) for g in gens(cox_ring(base))]
  if ("x" in gens_base_names) || ("y" in gens_base_names) || ("z" in gens_base_names)
    @vprint :WeierstrassModel 0 "Variable names duplicated between base and fiber coordinates.\n"
  end
  
  if completeness_check
    @req is_complete(base) "Base space must be complete"
  end
  
  ambient_space = _weierstrass_ambient_space_from_base(base)
  pw = _weierstrass_polynomial(f, g, cox_ring(ambient_space))
  model = WeierstrassModel(f, g, pw, toric_covered_scheme(base), toric_covered_scheme(ambient_space))
  set_attribute!(model, :base_fully_specified, true)
  return model
end


################################################
# 2: Constructors with toric scheme as base
################################################


@doc raw"""
    weierstrass_model(base::ToricCoveredScheme; completeness_check::Bool = true)

This method constructs a Weierstrass model over a given toric base
3-fold. The Weierstrass sections ``f`` and ``g`` are taken with (pseudo)random
coefficients.

# Examples
```jldoctest
julia> w = weierstrass_model(sample_toric_scheme(); completeness_check = false)
Weierstrass model over a concrete base
```
"""
weierstrass_model(base::ToricCoveredScheme; completeness_check::Bool = true) = weierstrass_model(underlying_toric_variety(base), completeness_check = completeness_check)


@doc raw"""
    weierstrass_model(f::MPolyRingElem, g::MPolyRingElem, base::ToricCoveredScheme; completeness_check::Bool = true)

This method operates analogously to `weierstrass_model(base::ToricCoveredScheme)`.
The only difference is that the Weierstrass sections ``f`` and ``g`` can be specified with non-generic values.

# Examples
```jldoctest
julia> base = sample_toric_scheme()
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[0, 1, 0], [0, 1, -1//2], [1, -1, -1], [0, 1, -1], [0, 0, 1], [0, 0, -1], [0, -1, 2], [0, -1, 1], [0, -1, 0], [0, -1, -1], [-1, 4, 0], [-1, 5, -1], [-1, 4, -1], [-1, 3, 1], [-1, 3, 0], [-1, 3, -1], [-1, 2, 2], [-1, 2, 1], [-1, 2, 0], [-1, 2, -1], [-1, 1, 3], [-1, 1, 2], [-1, 1, 1], [-1, 1, 0], [-1, 1, -1], [-1, 0, 4], [-1, 0, 3], [-1, 0, 2], [-1, 0, 1], [-1, 0, 0], [-1, 0, -1], [-1, -1, 5], [-1, -1, 4], [-1, -1, 3], [-1, -1, 2], [-1, -1, 1], [-1, -1, 0], [-1, -1, -1]]

julia> f = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(underlying_toric_variety(base))^4)]);

julia> g = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(underlying_toric_variety(base))^6)]);

julia> w = weierstrass_model(f, g, base; completeness_check = false)
Weierstrass model over a concrete base
```
"""
weierstrass_model(f::MPolyRingElem, g::MPolyRingElem, base::ToricCoveredScheme; completeness_check::Bool = true) = weierstrass_model(f, g, underlying_toric_variety(base), completeness_check = completeness_check)


################################################
# 3: Constructors with scheme as base
################################################

# Yet to come...
# This requires that the ai are stored as sections of the anticanonical bundle, and not "just" polynomials.
# -> Types to be generalized then.


################################################
# 4: Constructors without specified base
################################################

@doc raw"""
    weierstrass_model(weierstrass_f::MPolyRingElem, weierstrass_g::MPolyRingElem, auxiliary_base_ring::MPolyRing, d::Int)

This method constructs a Weierstrass model over a base space that is not
fully specified. The following example illustrates this approach.

# Examples
```jldoctest
julia> auxiliary_base_ring, (f, g, u) = QQ["f", "g", "u"];

julia> w = weierstrass_model(f, g, auxiliary_base_ring, 3)
Weierstrass model over a not fully specified base
```
"""
function weierstrass_model(weierstrass_f::MPolyRingElem, weierstrass_g::MPolyRingElem, auxiliary_base_ring::MPolyRing, d::Int)
  @req ((parent(weierstrass_f) == auxiliary_base_ring) && (parent(weierstrass_g) == auxiliary_base_ring)) "All Weierstrass sections must reside in the provided auxiliary base ring"
  @req d > 0 "The dimension of the base space must be positive"
  @req (ngens(auxiliary_base_ring) >= d) "We expect at least as many base variables as the desired base dimension"
  gens_base_names = [string(g) for g in gens(auxiliary_base_ring)]
  if ("x" in gens_base_names) || ("y" in gens_base_names) || ("z" in gens_base_names)
    @vprint :WeierstrassModel 0 "Variable names duplicated between base and fiber coordinates.\n"
  end
  
  # convert Weierstrass sections into polynomials of the auxiliary base
  auxiliary_base_space = _auxiliary_base_space(gens_base_names, d)
  S = cox_ring(auxiliary_base_space)
  ring_map = hom(auxiliary_base_ring, S, gens(S))
  f = ring_map(weierstrass_f)
  g = ring_map(weierstrass_g)
  
  # construct auxiliary ambient space
  auxiliary_ambient_space = _weierstrass_ambient_space_from_base(auxiliary_base_space)
  
  # compute model
  pw = _weierstrass_polynomial(f, g, cox_ring(auxiliary_ambient_space))
  model = WeierstrassModel(f, g, pw, toric_covered_scheme(auxiliary_base_space), toric_covered_scheme(auxiliary_ambient_space))
  set_attribute!(model, :base_fully_specified, false)
  return model
end


#######################################
# 5: Display
#######################################

function Base.show(io::IO, w::WeierstrassModel)
  if base_fully_specified(w)
    print(io, "Weierstrass model over a concrete base")
  else
    print(io, "Weierstrass model over a not fully specified base")
  end
end
