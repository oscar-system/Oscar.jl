################################################
# 1: The Julia type for GlobalWeierstrassModel
################################################

@attributes mutable struct GlobalWeierstrassModel
    poly_f::MPolyRingElem{QQFieldElem}
    poly_g::MPolyRingElem{QQFieldElem}
    pw::MPolyRingElem{QQFieldElem}
    toric_base_space::AbstractNormalToricVariety
    toric_ambient_space::AbstractNormalToricVariety
    Y4::ClosedSubvarietyOfToricVariety
    function GlobalWeierstrassModel(poly_f::MPolyRingElem{QQFieldElem},
                                    poly_g::MPolyRingElem{QQFieldElem},
                                    pw::MPolyRingElem{QQFieldElem},
                                    toric_base_space::AbstractNormalToricVariety,
                                    toric_ambient_space::AbstractNormalToricVariety,
                                    Y4::ClosedSubvarietyOfToricVariety)
        return new(poly_f, poly_g, pw, toric_base_space, toric_ambient_space, Y4)
    end
end


################################################
# 2: Constructors over specified bases
################################################

@doc Markdown.doc"""
    global_weierstrass_model(base::AbstractNormalToricVariety)

This method constructs a global Weierstrass model over a given toric base
3-fold. The Weierstrass sections ``f`` and ``g`` are taken with (pseudo)random
coefficients.

# Examples
```jldoctest
julia> w = global_weierstrass_model(test_base())
Global Weierstrass model over a concrete base

julia> dim(toric_ambient_space(w))
5
```
"""
function global_weierstrass_model(base::AbstractNormalToricVariety)
    toric_ambient_space = _ambient_space_from_base(base)
    (f,g) = _weierstrass_sections(base)
    pw = _weierstrass_polynomial(f, g, cox_ring(toric_ambient_space))
    Y4 = closed_subvariety_of_toric_variety(toric_ambient_space, [pw])
    model = GlobalWeierstrassModel(f, g, pw, base, toric_ambient_space, Y4)
    set_attribute!(model, :base_fully_specified, true)
    return model
end


@doc Markdown.doc"""
    global_weierstrass_model_over_projective_space()

This method constructs a global Weierstrass model over the 3-dimensional projective space.

# Examples
```jldoctest
julia> global_weierstrass_model_over_projective_space()
Global Weierstrass model over a concrete base
```
"""
global_weierstrass_model_over_projective_space() = global_weierstrass_model(projective_space(NormalToricVariety,3))


@doc Markdown.doc"""
    global_weierstrass_model(f::MPolyRingElem{QQFieldElem}, g::MPolyRingElem{QQFieldElem}, base::AbstractNormalToricVariety)

This method operates analogously to `global_weierstrass_model(base::AbstractNormalToricVariety)`.
The only difference is that the Weierstrass sections ``f`` and ``g`` can be specified with non-generic values.

# Examples
```jldoctest
julia> base = test_base()
Normal toric variety

julia> f = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^4)]);

julia> g = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)]);

julia> w = global_weierstrass_model(f, g, base)
Global Weierstrass model over a concrete base

julia> is_smooth(toric_ambient_space(w))
false
```
"""
function global_weierstrass_model(f::MPolyRingElem{QQFieldElem}, g::MPolyRingElem{QQFieldElem}, base::AbstractNormalToricVariety)
    @req ((parent(f) == cox_ring(base)) && (parent(g) == cox_ring(base))) "All Weierstrass sections must reside in the Cox ring of the base toric variety"
    toric_ambient_space = _ambient_space_from_base(base)
    pw = _weierstrass_polynomial(f, g, cox_ring(toric_ambient_space))
    Y4 = closed_subvariety_of_toric_variety(toric_ambient_space, [pw])
    model = GlobalWeierstrassModel(f, g, pw, base, toric_ambient_space, Y4)
    set_attribute!(model, :base_fully_specified, true)
    return model
end


################################################
# 3: Constructors over not fully specified bases
################################################

@doc Markdown.doc"""
    global_weierstrass_model(poly_f::MPolyRingElem{QQFieldElem}, poly_g::MPolyRingElem{QQFieldElem}, auxiliary_base_ring::MPolyRing, d::Int)

This method constructs a global Weierstrass model over a base space that is not
fully specified. The following example illustrates this approach.

# Examples
```jldoctest
julia> auxiliary_base_ring, (f, g, x) = QQ["f", "g", "x"];

julia> w = global_weierstrass_model(f, g, auxiliary_base_ring, 3)
Global Weierstrass model over a not fully specified base
```
"""
function global_weierstrass_model(poly_f::MPolyRingElem{QQFieldElem}, poly_g::MPolyRingElem{QQFieldElem}, auxiliary_base_ring::MPolyRing, d::Int)
    @req ((parent(poly_f) == auxiliary_base_ring) && (parent(poly_g) == auxiliary_base_ring)) "All Weierstrass sections must reside in the provided auxiliary base ring"
    @req d > 0 "The dimension of the base space must be positive"
    @req (ngens(auxiliary_base_ring) >= d) "We expect at least as many base variables as the desired base dimension"

    # convert Weierstrass sections into polynomials of the auxiliary base
    auxiliary_base_space = _auxiliary_base_space([string(k) for k in gens(auxiliary_base_ring)], d)
    S = cox_ring(auxiliary_base_space)
    ring_map = hom(auxiliary_base_ring, S, gens(S))
    f = ring_map(poly_f)
    g = ring_map(poly_g)

    # construct auxiliary ambient space
    auxiliary_ambient_space = _ambient_space_from_base(auxiliary_base_space)

    # compute model
    pw = _weierstrass_polynomial(f, g, cox_ring(auxiliary_ambient_space))
    Y4 = closed_subvariety_of_toric_variety(auxiliary_ambient_space, [pw])
    model = GlobalWeierstrassModel(f, g, pw, auxiliary_base_space, auxiliary_ambient_space, Y4)
    set_attribute!(model, :base_fully_specified, false)
    return model
end


#######################################
# 4: Display
#######################################

function Base.show(io::IO, w::GlobalWeierstrassModel)
    if base_fully_specified(w)
        print(io, "Global Weierstrass model over a concrete base")
    else
        print(io, "Global Weierstrass model over a not fully specified base")
    end
end
