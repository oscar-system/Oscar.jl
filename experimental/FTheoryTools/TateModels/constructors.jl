################################################
# 1: The Julia types for GlobalTateModel
################################################

@attributes mutable struct GlobalTateModel
    a1::MPolyRingElem{QQFieldElem}
    a2::MPolyRingElem{QQFieldElem}
    a3::MPolyRingElem{QQFieldElem}
    a4::MPolyRingElem{QQFieldElem}
    a6::MPolyRingElem{QQFieldElem}
    pt::MPolyRingElem{QQFieldElem}
    toric_base_space::AbstractNormalToricVariety
    toric_ambient_space::AbstractNormalToricVariety
    Y4::ClosedSubvarietyOfToricVariety
    function GlobalTateModel(a1::MPolyRingElem{QQFieldElem},
                            a2::MPolyRingElem{QQFieldElem},
                            a3::MPolyRingElem{QQFieldElem},
                            a4::MPolyRingElem{QQFieldElem},
                            a6::MPolyRingElem{QQFieldElem},
                            pt::MPolyRingElem{QQFieldElem},
                            toric_base_space::AbstractNormalToricVariety,
                            toric_ambient_space::AbstractNormalToricVariety,
                            Y4::ClosedSubvarietyOfToricVariety)
        return new(a1, a2, a3, a4, a6, pt, toric_base_space, toric_ambient_space, Y4)
    end
end


################################################
# 2: Constructors over specified bases
################################################

@doc Markdown.doc"""
    global_tate_model(base::AbstractNormalToricVariety)

This method constructs a global Tate model over a given toric base
3-fold. The Tate sections ``a_i`` are taken with (pseudo) random coefficients.

# Examples
```jldoctest
julia> t = global_tate_model(test_base())
Global Tate model over a concrete base

julia> is_smooth(toric_ambient_space(t))
false
```
"""
function global_tate_model(base::AbstractNormalToricVariety)
    toric_ambient_space = _ambient_space_from_base(base)
    (a1, a2, a3, a4, a6) = _tate_sections(base)
    pt = _tate_polynomial([a1, a2, a3, a4, a6], cox_ring(toric_ambient_space))
    Y4 = closed_subvariety_of_toric_variety(toric_ambient_space, [pt])
    model = GlobalTateModel(a1, a2, a3, a4, a6, pt, base, toric_ambient_space, Y4)
    set_attribute!(model, :base_fully_specified, true)
    return model
end


@doc Markdown.doc"""
    global_tate_model_over_projective_space()

This method constructs a global Tate model over the 3-dimensional projective space.

# Examples
```jldoctest
julia> global_tate_model_over_projective_space()
Global Tate model over a concrete base
```
"""
global_tate_model_over_projective_space() = global_tate_model(projective_space(NormalToricVariety,3))


@doc Markdown.doc"""
    global_tate_model(ais::Vector{T}, base::AbstractNormalToricVariety) where {T<:MPolyRingElem{QQFieldElem}}

This method operates analogously to `global_tate_model(base::AbstractNormalToricVariety)`.
The only difference is that the Tate sections ``a_i`` can be specified with non-generic values.

# Examples
```jldoctest
julia> base = test_base()
Normal toric variety

julia> a1 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base))]);

julia> a2 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^2)]);

julia> a3 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^3)]);

julia> a4 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^4)]);

julia> a6 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)]);

julia> t = global_tate_model([a1, a2, a3, a4, a6], base)
Global Tate model over a concrete base

julia> is_smooth(toric_ambient_space(t))
false
```
"""
function global_tate_model(ais::Vector{T}, base::AbstractNormalToricVariety) where {T<:MPolyRingElem{QQFieldElem}}
    if length(ais) != 5
        throw(ArgumentError("We require exactly 5 Tate sections"))
    end
    if any(k -> parent(k) != cox_ring(base), ais)
        throw(ArgumentError("All Tate sections must reside in the Cox ring of the base toric variety"))
    end
    toric_ambient_space = _ambient_space_from_base(base)
    pt = _tate_polynomial(ais, cox_ring(toric_ambient_space))
    Y4 = closed_subvariety_of_toric_variety(toric_ambient_space, [pt])
    model = GlobalTateModel(ais[1], ais[2], ais[3], ais[4], ais[5], pt, base, toric_ambient_space, Y4)
    set_attribute!(model, :base_fully_specified, true)
    return model
end


################################################
# 3: Constructors over not fully specified bases
################################################

@doc Markdown.doc"""
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

julia> tate_polynomial(t)
-a10*x*y*z + a21*w*x^2*z^2 - a32*w^2*y*z^3 + a43*w^3*x*z^4 + a65*w^5*z^6 + x^3 - y^2

julia> toric_base_space(t)
[ Info: Base space was not fully specified. Returning AUXILIARY base space.
Normal toric variety

julia> toric_ambient_space(t)
[ Info: Base space was not fully specified. Returning AUXILIARY ambient space.
Normal, simplicial toric variety

julia> dim(toric_ambient_space(t))
5
```
"""
function global_tate_model(ais::Vector{T}, auxiliary_base_ring::MPolyRing, d::Int) where {T<:MPolyRingElem{QQFieldElem}}
    if length(ais) != 5
        throw(ArgumentError("We expect exactly 5 Tate sections"))
    end
    if any(k -> parent(k) != auxiliary_base_ring, ais)
        throw(ArgumentError("All Tate sections must reside in the provided auxiliary base ring"))
    end
    if d <= 0
        throw(ArgumentError("The dimension of the base space must be positive"))
    end
    if ngens(auxiliary_base_ring) < d
        throw(ArgumentError("We expect at least as many base variables as the desired base dimension"))
    end

    # convert Tate sections into polynomials of the auxiliary base
    auxiliary_base_space = _auxiliary_base_space([string(k) for k in gens(auxiliary_base_ring)], d)
    S = cox_ring(auxiliary_base_space)
    ring_map = hom(auxiliary_base_ring, S, gens(S))
    (a1, a2, a3, a4, a6) = [ring_map(k) for k in ais]

    # construct model
    auxiliary_ambient_space = _ambient_space_from_base(auxiliary_base_space)
    pt = _tate_polynomial([a1, a2, a3, a4, a6], cox_ring(auxiliary_ambient_space))
    Y4 = closed_subvariety_of_toric_variety(auxiliary_ambient_space, [pt])
    model = GlobalTateModel(a1, a2, a3, a4, a6, pt, auxiliary_base_space, auxiliary_ambient_space, Y4)
    set_attribute!(model, :base_fully_specified, false)
    return model
end


################################################
# 4: Display
################################################

function Base.show(io::IO, t::GlobalTateModel)
    if base_fully_specified(t)
        join(io, "Global Tate model over a concrete base")
    else
        join(io, "Global Tate model over a not fully specified base")
    end
end
