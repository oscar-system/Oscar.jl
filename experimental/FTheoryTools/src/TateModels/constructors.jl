################################################
# 1: Constructors with toric variety as base
################################################

@doc raw"""
    global_tate_model(base::AbstractNormalToricVariety; completeness_check::Bool = true)

This method constructs a global Tate model over a given toric base
3-fold. The Tate sections ``a_i`` are taken with (pseudo) random coefficients.

# Examples
```jldoctest
julia> t = global_tate_model(sample_toric_variety(); completeness_check = false)
Global Tate model over a concrete base
```
"""
global_tate_model(base::AbstractNormalToricVariety; completeness_check::Bool = true) = global_tate_model(_tate_sections(base), base; completeness_check = completeness_check)


@doc raw"""
    global_tate_model(ais::Vector{T}, base::AbstractNormalToricVariety; completeness_check::Bool = true) where {T<:MPolyRingElem{QQFieldElem}}

This method operates analogously to `global_tate_model(base::AbstractNormalToricVariety)`.
The only difference is that the Tate sections ``a_i`` can be specified with non-generic values.

# Examples
```jldoctest
julia> base = sample_toric_variety()
Normal toric variety

julia> a1 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base))]);

julia> a2 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^2)]);

julia> a3 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^3)]);

julia> a4 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^4)]);

julia> a6 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)]);

julia> t = global_tate_model([a1, a2, a3, a4, a6], base; completeness_check = false)
Global Tate model over a concrete base
```
"""
function global_tate_model(ais::Vector{T}, base::AbstractNormalToricVariety; completeness_check::Bool = true) where {T<:MPolyRingElem{QQFieldElem}}
  if completeness_check
    @req is_complete(base) "Base space must be complete"
  end
  @req length(ais) == 5 "We require exactly 5 Tate sections"
  @req all(k -> parent(k) == cox_ring(base), ais) "All Tate sections must reside in the Cox ring of the base toric variety"
  ambient_space = _ambient_space_from_base(base)
  pt = _tate_polynomial(ais, cox_ring(ambient_space))
  model = GlobalTateModel(ais[1], ais[2], ais[3], ais[4], ais[5], pt, toric_covered_scheme(base), toric_covered_scheme(ambient_space))
  set_attribute!(model, :base_fully_specified, true)
  return model
end


################################################
# 2: Constructors with toric scheme as base
################################################


@doc raw"""
    global_tate_model(base::ToricCoveredScheme; completeness_check::Bool = true)

This method constructs a global Tate model over a given toric scheme base
3-fold. The Tate sections ``a_i`` are taken with (pseudo) random coefficients.

# Examples
```jldoctest
julia> t = global_tate_model(sample_toric_scheme(); completeness_check = false)
Global Tate model over a concrete base
```
"""
global_tate_model(base::ToricCoveredScheme; completeness_check::Bool = true) = global_tate_model(underlying_toric_variety(base), completeness_check = completeness_check)


@doc raw"""
    global_tate_model(ais::Vector{T}, base::ToricCoveredScheme; completeness_check::Bool = true) where {T<:MPolyRingElem{QQFieldElem}}

This method operates analogously to `global_tate_model(base::ToricCoveredScheme)`.
The only difference is that the Tate sections ``a_i`` can be specified with non-generic values.

# Examples
```jldoctest
julia> base = sample_toric_scheme()
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[0, 1, 0], [0, 1, -1//2], [1, -1, -1], [0, 1, -1], [0, 0, 1], [0, 0, -1], [0, -1, 2], [0, -1, 1], [0, -1, 0], [0, -1, -1], [-1, 4, 0], [-1, 5, -1], [-1, 4, -1], [-1, 3, 1], [-1, 3, 0], [-1, 3, -1], [-1, 2, 2], [-1, 2, 1], [-1, 2, 0], [-1, 2, -1], [-1, 1, 3], [-1, 1, 2], [-1, 1, 1], [-1, 1, 0], [-1, 1, -1], [-1, 0, 4], [-1, 0, 3], [-1, 0, 2], [-1, 0, 1], [-1, 0, 0], [-1, 0, -1], [-1, -1, 5], [-1, -1, 4], [-1, -1, 3], [-1, -1, 2], [-1, -1, 1], [-1, -1, 0], [-1, -1, -1]]

julia> a1 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(underlying_toric_variety(base)))]);

julia> a2 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(underlying_toric_variety(base))^2)]);

julia> a3 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(underlying_toric_variety(base))^3)]);

julia> a4 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(underlying_toric_variety(base))^4)]);

julia> a6 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(underlying_toric_variety(base))^6)]);

julia> t = global_tate_model([a1, a2, a3, a4, a6], base; completeness_check = false)
Global Tate model over a concrete base
```
"""
global_tate_model(ais::Vector{T}, base::ToricCoveredScheme; completeness_check::Bool = true) where {T<:MPolyRingElem{QQFieldElem}} = global_tate_model(ais, underlying_toric_variety(base); completeness_check = completeness_check)


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
    global_tate_model(ais::Vector{T}, auxiliary_base_ring::MPolyRing, d::Int) where {T<:MPolyRingElem{QQFieldElem}}

This method constructs a global Tate model over a base space that is not
fully specified. The following example exemplifies this approach.

# Examples
```jldoctest
julia> auxiliary_base_ring, (a10, a21, a32, a43, a65, w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];

julia> a1 = a10;

julia> a2 = a21 * w;

julia> a3 = a32 * w^2;

julia> a4 = a43 * w^3;

julia> a6 = a65 * w^5;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = global_tate_model(ais, auxiliary_base_ring, 3)
Global Tate model over a not fully specified base
```
"""
function global_tate_model(ais::Vector{T}, auxiliary_base_ring::MPolyRing, d::Int) where {T<:MPolyRingElem{QQFieldElem}}
  @req length(ais) == 5 "We expect exactly 5 Tate sections"
  @req all(k -> parent(k) == auxiliary_base_ring, ais) "All Tate sections must reside in the provided auxiliary base ring"
  @req d > 0 "The dimension of the base space must be positive"
  @req ngens(auxiliary_base_ring) >= d "We expect at least as many base variables as the desired base dimension"
  
  # convert Tate sections into polynomials of the auxiliary base
  auxiliary_base_space = _auxiliary_base_space([string(k) for k in gens(auxiliary_base_ring)], d)
  S = cox_ring(auxiliary_base_space)
  ring_map = hom(auxiliary_base_ring, S, gens(S))
  (a1, a2, a3, a4, a6) = [ring_map(k) for k in ais]
  
  # construct model
  auxiliary_ambient_space = _ambient_space_from_base(auxiliary_base_space)
  pt = _tate_polynomial([a1, a2, a3, a4, a6], cox_ring(auxiliary_ambient_space))
  model = GlobalTateModel(a1, a2, a3, a4, a6, pt, toric_covered_scheme(auxiliary_base_space), toric_covered_scheme(auxiliary_ambient_space))
  set_attribute!(model, :base_fully_specified, false)
  return model
end


################################################
# 5: Display
################################################

function Base.show(io::IO, t::GlobalTateModel)
  if base_fully_specified(t)
    print(io, "Global Tate model over a concrete base")
  else
    print(io, "Global Tate model over a not fully specified base")
  end
end
