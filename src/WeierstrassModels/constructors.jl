################################################
# 1: The Julia type for GlobalWeierstrassModel
################################################

@attributes mutable struct GlobalWeierstrassModel
    poly_f::MPolyElem{fmpq}
    poly_g::MPolyElem{fmpq}
    pw::MPolyElem{fmpq}
    toric_base_space::Oscar.AbstractNormalToricVariety
    toric_ambient_space::Oscar.AbstractNormalToricVariety
    Y4::Oscar.ClosedSubvarietyOfToricVariety
    function GlobalWeierstrassModel(poly_f::MPolyElem{fmpq},
                                    poly_g::MPolyElem{fmpq},
                                    pw::MPolyElem{fmpq},
                                    toric_base_space::Oscar.AbstractNormalToricVariety,
                                    toric_ambient_space::Oscar.AbstractNormalToricVariety,
                                    Y4::Oscar.ClosedSubvarietyOfToricVariety)
        return new(poly_f, poly_g, pw, toric_base_space, toric_ambient_space, Y4)
    end
end
export GlobalWeierstrassModel


################################################
# 2: Constructors over specified bases
################################################

@doc Markdown.doc"""
    GlobalWeierstrassModel(base::Oscar.AbstractNormalToricVariety)

This method constructs a global Weierstrass model over a given toric base
3-fold. The Weierstrass sections ``f`` and ``g`` are taken with (pseudo)
random coefficients.

# Examples
```jldoctest
julia> using Oscar

julia> w = GlobalWeierstrassModel(TestBase())
A global Weierstrass model over a concrete base

julia> is_smooth(toric_ambient_space(w))
false
```
"""
function GlobalWeierstrassModel(base::Oscar.AbstractNormalToricVariety)
    toric_ambient_space = _ambient_space_from_base(base)
    (f,g) = _weierstrass_sections(base)
    pw = _weierstrass_polynomial(f, g, cox_ring(toric_ambient_space))
    Y4 = Oscar.ClosedSubvarietyOfToricVariety(toric_ambient_space, [pw])
    model = GlobalWeierstrassModel(f, g, pw, base, toric_ambient_space, Y4)
    set_attribute!(model, :base_fully_specified, true)
    return model
end
export GlobalWeierstrassModel


@doc Markdown.doc"""
    GlobalWeierstrassModelOverProjectiveSpace()

This method constructs a global Weierstrass model over the 3-dimensional projective space.

# Examples
```jldoctest
julia> using Oscar

julia> GlobalWeierstrassModelOverProjectiveSpace()
A global Weierstrass model over a concrete base
```
"""
GlobalWeierstrassModelOverProjectiveSpace() = GlobalWeierstrassModel(projective_space(NormalToricVariety,3))
export GlobalWeierstrassModelOverProjectiveSpace



@doc Markdown.doc"""
    GlobalWeierstrassModel(f::MPolyElem{fmpq}, g::MPolyElem{fmpq}, base::Oscar.AbstractNormalToricVariety)

This method operates analogously to `GenericGlobalWeierstrassModel(base::Oscar.AbstractNormalToricVariety)`.
The only difference is that the Weierstrass sections ``f`` and ``g`` can be specified with non-generic values.

# Examples
```jldoctest
julia> using Oscar

julia> base = TestBase()
A normal toric variety

julia> f = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^4)]);

julia> g = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)]);

julia> w = GlobalWeierstrassModel(f, g, base)
A global Weierstrass model over a concrete base

julia> is_smooth(toric_ambient_space(w))
false
```
"""
function GlobalWeierstrassModel(f::MPolyElem{fmpq}, g::MPolyElem{fmpq}, base::Oscar.AbstractNormalToricVariety)
    if (parent(f) != cox_ring(base)) || (parent(g) != cox_ring(base))
        throw(ArgumentError("All Weierstrass sections must reside in the Cox ring of the base toric variety"))
    end
    toric_ambient_space = _ambient_space_from_base(base)
    pw = _weierstrass_polynomial(f, g, cox_ring(toric_ambient_space))
    Y4 = Oscar.ClosedSubvarietyOfToricVariety(toric_ambient_space, [pw])
    model = GlobalWeierstrassModel(f, g, pw, base, toric_ambient_space, Y4)
    set_attribute!(model, :base_fully_specified, true)
    return model
end
export GlobalWeierstrassModel


################################################
# 3: Constructors over not fully specified bases
################################################

@doc Markdown.doc"""
    GlobalWeierstrassModel(poly_f::MPolyElem{fmpq}, poly_g::MPolyElem{fmpq}, auxiliary_base_ring::MPolyRing, d::Int)

This method constructs a global Weierstrass model over a base space that is not
fully specified. The following example illustrates this approach.

# Examples
```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (f, g, x) = QQ["f", "g", "x"];

julia> w = GlobalWeierstrassModel(f, g, auxiliary_base_ring, 3)
A global Weierstrass model over a not fully specified base

julia> weierstrass_polynomial(w)
f*x*z^4 + g*z^6 + x^3 - y^2

julia> toric_base_space(w)
[ Info: Base space was not fully specified. Returning AUXILIARY base space.
A normal toric variety

julia> toric_ambient_space(w)
[ Info: Base space was not fully specified. Returning AUXILIARY ambient space.
A normal, simplicial toric variety

julia> dim(toric_ambient_space(w))
5
```
"""
function GlobalWeierstrassModel(poly_f::MPolyElem{fmpq}, poly_g::MPolyElem{fmpq}, auxiliary_base_ring::MPolyRing, d::Int)
    if (parent(poly_f) != auxiliary_base_ring) || (parent(poly_g) != auxiliary_base_ring)
        throw(ArgumentError("All Weierstrass sections must reside in the provided auxiliary base ring"))
    end
    if d <= 0
        throw(ArgumentError("The dimension of the base space must be positive"))
    end
    if length([string(k) for k in gens(auxiliary_base_ring)]) < d
        throw(ArgumentError("We expect at least as many base variables as the desired base dimension"))
    end

    # convert Tate sections into polynomials of the auxiliary base
    auxiliary_base_space = _auxiliary_base_space([string(k) for k in gens(auxiliary_base_ring)], d)
    S = cox_ring(auxiliary_base_space)
    ring_map = hom(auxiliary_base_ring, S, [gens(S)[i] for i in 1:length(gens(S))])
    f = ring_map(poly_f)
    g = ring_map(poly_g)

    # construct auxiliary ambient space
    auxiliary_ambient_space = _ambient_space_from_base(auxiliary_base_space)

    # compute model
    pw = _weierstrass_polynomial(f, g, cox_ring(auxiliary_ambient_space))
    Y4 = ClosedSubvarietyOfToricVariety(auxiliary_ambient_space, [pw])
    model = GlobalWeierstrassModel(f, g, pw, auxiliary_base_space, auxiliary_ambient_space, Y4)
    set_attribute!(model, :base_fully_specified, false)
    return model
end
export GlobalWeierstrassModel


#######################################
# 4: Display
#######################################

function Base.show(io::IO, w::GlobalWeierstrassModel)
    if base_fully_specified(w)
        join(io, "A global Weierstrass model over a concrete base")
    else
        join(io, "A global Weierstrass model over a not fully specified base")
    end
end
