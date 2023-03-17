################################################################
# 1: Construct auxiliary base space
################################################################

function _auxiliary_base_space(variable_names::Vector{String}, d::Int)
    ray_gens = [[if i==j 1 else 0 end for j in 1:length(variable_names)] for i in 1:length(variable_names)]
    max_cones = Oscar.Hecke.subsets(length(variable_names), d)
    auxiliary_base_space = NormalToricVariety(ray_gens, max_cones)
    set_coordinate_names(auxiliary_base_space, variable_names)
    return auxiliary_base_space
end

################################################################
# 2: Construct ambient space from given base
################################################################

function _ambient_space_from_base(base::Oscar.AbstractNormalToricVariety)

    # Extract information about the toric base
    base_rays = matrix(ZZ, rays(base))
    base_cones = matrix(ZZ, ray_indices(maximal_cones(base)))

    # Construct the rays for the fibre ambient space
    xray = [0 for i in 1:ncols(base_rays)+2]
    yray = [0 for i in 1:ncols(base_rays)+2]
    zray = [0 for i in 1:ncols(base_rays)+2]
    xray[ncols(base_rays)+1] = 1
    yray[ncols(base_rays)+2] = 1
    zray[ncols(base_rays)+1] = -2
    zray[ncols(base_rays)+2] = -3

    # Construct the rays of the toric ambient space
    ambient_space_rays = hcat([r for r in base_rays], [-2 for i in 1:nrows(base_rays)], [-3 for i in 1:nrows(base_rays)])
    ambient_space_rays = vcat(ambient_space_rays, transpose(xray), transpose(yray), transpose(zray))

    # Construct the incidence matrix for the maximal cones of the ambient space
    ambient_space_max_cones = []
    for i in 1:nrows(base_cones)
        push!(ambient_space_max_cones, [k for k in hcat([b for b in base_cones[i,:]], [1 1 0])])
        push!(ambient_space_max_cones, [k for k in hcat([b for b in base_cones[i,:]], [1 0 1])])
        push!(ambient_space_max_cones, [k for k in hcat([b for b in base_cones[i,:]], [0 1 1])])
    end
    ambient_space_max_cones = IncidenceMatrix(vcat(ambient_space_max_cones...))

    # Construct and return the ambient space
    toric_ambient_space = NormalToricVariety(PolyhedralFan(ambient_space_rays, ambient_space_max_cones; non_redundant = true))
    set_coordinate_names(toric_ambient_space, vcat([string(k) for k in gens(cox_ring(base))], ["x", "y", "z"]))
    return toric_ambient_space

end


################################################################
# 3: Construct the Weierstrass polynomial
################################################################

function _weierstrass_sections(base::Oscar.AbstractNormalToricVariety)
    f = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^4)])
    g = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)])
    return [f, g]
end

function _weierstrass_polynomial(base::Oscar.AbstractNormalToricVariety, S::MPolyRing_dec{fmpq, FmpqMPolyRing})
    (f, g) = _weierstrass_sections(base)
    return _weierstrass_polynomial(S, f, g)
end

function _weierstrass_polynomial(f::MPolyElem{fmpq}, g::MPolyElem{fmpq}, S::MPolyRing_dec{fmpq, FmpqMPolyRing})
    x = gens(S)[length(gens(S))-2]
    y = gens(S)[length(gens(S))-1]
    z = gens(S)[length(gens(S))]
    ring_map = hom(parent(f), S, [gens(S)[i] for i in 1:length(gens(S))-3])
    return x^3 - y^2 + ring_map(f)*x*z^4 + ring_map(g)*z^6
end


################################################################
# 4: Construct the Tate polynomial
################################################################

function _tate_sections(base::Oscar.AbstractNormalToricVariety)
    a1 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base))])
    a2 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^2)])
    a3 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^3)])
    a4 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^4)])
    a6 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)])
    return [a1, a2, a3, a4, a6]
end

function _tate_polynomial(base::Oscar.AbstractNormalToricVariety, S::MPolyRing_dec{fmpq, FmpqMPolyRing})
    (a1, a2, a3, a4, a6) = _tate_sections(base)
    return _tate_polynomial(S, [a1, a2, a3, a4, a6])
end

function _tate_polynomial(ais::Vector{<:MPolyElem{fmpq}}, S::MPolyRing_dec{fmpq, FmpqMPolyRing})
    x = gens(S)[length(gens(S))-2]
    y = gens(S)[length(gens(S))-1]
    z = gens(S)[length(gens(S))]
    ring_map = hom(parent(ais[1]), S, [gens(S)[i] for i in 1:length(gens(S))-3])
    (a1, a2, a3, a4, a6) = [ring_map(k) for k in ais]
    return x^3 - y^2 - x*y*z*a1 + x^2*z^2*a2 - y*z^3*a3 + x*z^4*a4 + z^6*a6
end


################################################################
# 5: A base space for efficient testing
################################################################

@doc Markdown.doc"""
    test_base()

This method constructs a 3-dimensional toric variety, which we
use for efficient testing of the provided functionality.
"""
function test_base()
    rays = [-1 -1 -1; -1 -1 0; -1 -1 1; -1 -1 2; -1 -1 3; -1 -1 4; -1 -1 5; -1 0 -1; -1 0 0; -1 0 1; -1 0 2; -1 0 3; -1 0 4; -1 1 -1; -1 1 0; -1 1 1; -1 1 2; -1 1 3; -1 2 -1; -1 2 0; -1 2 1; -1 2 2; -1 3 -1; -1 3 0; -1 3 1; -1 4 -1; -1 4 0; -1 5 -1; 0 -1 -1; 0 -1 0; 0 -1 1; 0 -1 2; 0 0 -1; 0 0 1; 0 1 -1; 0 1 0; 0 2 -1; 1 -1 -1]
    cones = IncidenceMatrix([[36, 37, 38], [35, 37, 38], [34, 36, 38], [33, 35, 38], [32, 34, 38], [31, 32, 38], [30, 31, 38], [29, 33, 38], [29, 30, 38], [27, 28, 37], [26, 28, 37], [26, 27, 28], [25, 36, 37], [25, 27, 37], [24, 26, 27], [24, 25, 27], [23, 26, 37], [23, 24, 26], [22, 34, 36], [22, 25, 36], [21, 24, 25], [21, 22, 25], [20, 23, 24], [20, 21, 24], [19, 35, 37], [19, 23, 37], [19, 20, 23], [18, 32, 34], [18, 22, 34], [17, 21, 22], [17, 18, 22], [16, 20, 21], [16, 17, 21], [15, 19, 20], [15, 16, 20], [14, 33, 35], [14, 19, 35], [14, 15, 19], [13, 18, 32], [12, 17, 18], [12, 13, 18], [11, 16, 17], [11, 12, 17], [10, 15, 16], [10, 11, 16], [9, 14, 15], [9, 10, 15], [8, 29, 33], [8, 14, 33], [8, 9, 14], [7, 13, 32], [6, 12, 13], [6, 7, 32], [6, 7, 13], [5, 11, 12], [5, 6, 32], [5, 6, 12], [4, 31, 32], [4, 10, 11], [4, 5, 32], [4, 5, 11], [3, 30, 31], [3, 9, 10], [3, 4, 31], [3, 4, 10], [2, 29, 30], [2, 8, 9], [2, 3, 30], [2, 3, 9], [1, 8, 29], [1, 2, 29], [1, 2, 8]])
    return NormalToricVariety(PolyhedralFan(rays, cones))
end
export test_base


################################################################
# 6: Check if an ideal/subvariety is nontrivial
################################################################

_is_nontrivial(id::MPolyIdeal{T}, irr::MPolyIdeal{T}) where {T<:MPolyElem{fmpq}} = !is_one(id) && !is_one(saturation(id, irr))
export _is_nontrivial


################################################################
# 7: Compute singularity Kodaira type and refined Tate type
################################################################

_count_factors(poly::fmpq_mpoly) = mapreduce(p -> p[end], +, absolute_primary_decomposition(ideal([poly])))

function _kodaira_type(id::MPolyIdeal{T}, f::T, g::T, d::T, ords::Tuple{Int64, Int64, Int64}) where {T<:MPolyElem_dec{fmpq, fmpq_mpoly}}
    f_ord = ords[1]
    g_ord = ords[2]
    d_ord = ords[3]

    if d_ord == 0
        kod_type = "I_0"
    elseif d_ord == 1 && f_ord == 0 && g_ord == 0
        kod_type = "I_1"
    elseif d_ord == 2 && g_ord == 1 && f_ord >= 1
        kod_type = "II"
    elseif d_ord == 3 && f_ord == 1 && g_ord >= 2
        kod_type = "III"
    elseif d_ord == 9 && f_ord == 3 && g_ord >= 5
        kod_type = "III^*"
    elseif d_ord == 10 && g_ord == 5 && f_ord >= 4
        kod_type = "II^*"
    elseif d_ord >= 12 && f_ord >= 4 && g_ord >= 6
        kod_type = "Non-minimal"
    else
        R = parent(f)
        S, (_psi, ) = PolynomialRing(QQ, ["_psi"; [string(v) for v in gens(R)]], cached = false)
        ring_map = hom(R, S, gens(S)[2:end])
        poly_f = ring_map(f)
        poly_g = ring_map(g)
        poly_d = ring_map(d)
        locus = ring_map(gens(id)[1])

        if f_ord == 0 && g_ord == 0
            monodromy_poly = _psi^2 + divexact(evaluate(9 * poly_g, [locus], [0]), evaluate(2 * poly_f, [locus], [0]))
            if _count_factors(monodromy_poly) == 2
                kod_type = "Split I_$d_ord"
            else
                kod_type = "Non-split I_$d_ord"
            end
        elseif d_ord == 4 && g_ord == 2 && f_ord >= 2
            monodromy_poly = _psi^2 - evaluate(divexact(poly_g, locus^2), [locus], [0])
            if _count_factors(monodromy_poly) == 2
                kod_type = "Split IV"
            else
                kod_type = "Non-split IV"
            end
        elseif d_ord == 6 && f_ord >= 2 && g_ord >= 3
            monodromy_poly =  _psi^3 + _psi * evaluate(divexact(poly_f, locus^2), [locus], [0]) + evaluate(divexact(poly_g, locus^3), [locus], [0])
            num_facs = _count_factors(monodromy_poly)
            if num_facs == 3
                kod_type = "Split I^*_0"
            elseif num_facs == 2
                kod_type = "Semi-split I^*_0"
            else
                kod_type = "Non-split I^*_0"
            end
        elseif f_ord == 2 && g_ord == 3 && d_ord >= 7 && d_ord % 2 == 1
            monodromy_poly = _psi^2 + divexact(evaluate(divexact(poly_d, locus^d_ord) * divexact(2 * poly_f, locus^2)^3, [locus], [0]), 4 * evaluate(divexact(9 * poly_g, locus^3), [locus], [0])^3)
            if _count_factors(monodromy_poly) == 2
                kod_type = "Split I^*_$(d_ord - 6)"
            else
                kod_type = "Non-split I^*_$(d_ord - 6)"
            end
        elseif f_ord == 2 && g_ord == 3 && d_ord >= 8 && d_ord % 2 == 0
            monodromy_poly = _psi^2 + divexact(evaluate(divexact(poly_d, locus^d_ord) * divexact(2 * poly_f, locus^2)^2, [locus], [0]), evaluate(divexact(9 * poly_g, locus^3), [locus], [0])^2)
            if _count_factors(monodromy_poly) == 2
                kod_type = "Split I^*_$(d_ord - 6)"
            else
                kod_type = "Non-split I^*_$(d_ord - 6)"
            end
        elseif d_ord == 8 && g_ord == 4 && f_ord >= 3
            monodromy_poly = _psi^2 - evaluate(divexact(poly_g, locus^4), [locus], [0])
            if _count_factors(monodromy_poly) == 2
                kod_type = "Split IV^*"
            else
                kod_type = "Non-split IV^*"
            end
        else
            kod_type = "Unrecognized"
        end
    end
    
    return kod_type
end


################################################################
# 8: Blowups
################################################################

function _blowup_global(id::MPolyIdeal{fmpq_mpoly}, center::MPolyIdeal{fmpq_mpoly}, irr::MPolyIdeal{fmpq_mpoly}, sri::MPolyIdeal{fmpq_mpoly}, lin::MPolyIdeal{fmpq_mpoly}; index::Integer = 1)
    # @warn "The function _blowup_global is experimental; absence of bugs and proper results are not guaranteed"

    R = base_ring(id)
    center_size = length(gens(center))

    # Various sanity checks
    if is_zero(center)
        throw(ArgumentError("The blowup center must be non-empty"))
    end
    # if !is_subset(id, center)
    #     throw(ArgumentError("The ideal of the blowup center must contain the ideal to be blown up"))
    # end
    if base_ring(irr) != R
        throw(ArgumentError("The given irrelevant ideal must share the base ring of the ideal to be blown up"))
    end
    if base_ring(sri) != R
        throw(ArgumentError("The given Stanley–Reisner ideal must share the base ring of the ideal to be blown up"))
    end
    if length(gens(base_ring(lin))) != length(gens(R))
        throw(ArgumentError("The base ring of ideal of linear relations must have the same number of generators as the base ring of the ideal to be blown up"))
    end

    # Make sure the ideal of linear relations has the same base ring as the others
    lin = ideal(map(hom(base_ring(lin), R, collect(1:length(gens(R)))), gens(lin)))

    # Create new base ring for the blown up ideal and a map between the rings
    S, S_gens = PolynomialRing(QQ, [string("e_", index); [string("b_", index, "_", i) for i in 1:center_size]; [string(v) for v in gens(R)]], cached = false)
    (_e, new_coords...) = S_gens[1:center_size + 1]
    ring_map = hom(R, S, S_gens[center_size + 2:end])

    # Compute the total transform
    center_gens_S = map(ring_map, gens(center))
    total_transform = ideal(map(ring_map, gens(id))) + ideal([new_coords[i] * _e - center_gens_S[i] for i in 1:center_size])

    # Compute the exceptional locus and strict transform, checking for crepancy
    # Could alternatively replace _e with center_gens_S in the exceptional locus here, then take the
    #   primary decomposition and remove parts whose saturation by the irrelevant ideal is the whole ring
    exceptional_ideal = total_transform + ideal([_e])
    strict_transform, exceptional_factor = saturation_with_index(total_transform, exceptional_ideal)
    crepant = (exceptional_factor == center_size - 1)

    # Compute the new irrelevant ideal, SRI, and ideal of linear relations
    # These may need to be changed after reintroducing e
    new_irr = ideal(map(ring_map, gens(irr))) * ideal(new_coords)
    new_sri = ideal(map(ring_map, gens(sri))) + ideal([prod(new_coords)])
    new_lin = ideal(map(ring_map, gens(lin))) + ideal([g - new_coords[end] for g in new_coords[1:end - 1]])

    return total_transform, strict_transform, exceptional_ideal, crepant, new_irr, new_sri, new_lin, S, S_gens, ring_map
end
_blowup_global(id::T, center::T, irr::T, sri::T, lin::MPolyIdeal{fmpq_mpoly}; index::Integer = 1) where {T<:MPolyIdeal{MPolyElem_dec{fmpq, fmpq_mpoly}}} = _blowup_global(ideal(map(g -> g.f, gens(id))), ideal(map(g -> g.f, gens(center))), ideal(map(g -> g.f, gens(irr))), ideal(map(g -> g.f, gens(sri))), lin, index = index)
export _blowup_global

function _blowup_global_sequence(id::MPolyIdeal{fmpq_mpoly}, centers::Vector{<:Vector{<:Integer}}, irr::MPolyIdeal{fmpq_mpoly}, sri::MPolyIdeal{fmpq_mpoly}, lin::MPolyIdeal{fmpq_mpoly}; index::Integer = 1)
    # @warn "The function _blowup_global_sequence is experimental; absence of bugs and proper results are not guaranteed"

    (cur_strict_transform, cur_irr, cur_sri, cur_lin, cur_S, cur_S_gens, cur_index) = (id, irr, sri, lin, base_ring(id), gens(base_ring((id))), index)
    crepant = true
    ring_map = hom(cur_S, cur_S, cur_S_gens) # Identity map

    exceptionals = MPolyIdeal{<:MPolyElem{fmpq}}[]
    for center in centers
        if !all(ind -> 1 <= ind <= length(cur_S_gens), center)
            throw(ArgumentError("The given indices for the center generators are out of bounds"))
        end

        (_, cur_strict_transform, cur_ex, cur_crep, cur_irr, cur_sri, cur_lin, cur_S, cur_S_gens, cur_ring_map) = _blowup_global(cur_strict_transform, ideal(map(ind -> cur_S_gens[ind], center)), cur_irr, cur_sri, cur_lin, index = cur_index)

        map!(cur_ring_map, exceptionals, exceptionals)
        push!(exceptionals, cur_ex)

        crepant = crepant && cur_crep 

        ring_map = compose(ring_map, cur_ring_map)

        cur_index += 1
    end

    return cur_strict_transform, exceptionals, crepant, cur_irr, cur_sri, cur_lin, cur_S, cur_S_gens, ring_map
end
_blowup_global_sequence(id::T, centers::Vector{<:Vector{<:Integer}}, irr::T, sri::T, lin::MPolyIdeal{fmpq_mpoly}; index::Integer = 1) where {T<:MPolyIdeal{MPolyElem_dec{fmpq, fmpq_mpoly}}} = _blowup_global_sequence(ideal(map(g -> g.f, gens(id))), centers, ideal(map(g -> g.f, gens(irr))), ideal(map(g -> g.f, gens(sri))), lin, index = index)
export _blowup_global_sequence
