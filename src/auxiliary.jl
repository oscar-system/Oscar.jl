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
    TestBase()

This method constructs a 3-dimensional toric variety, which we
use for efficient testing of the provided functionality.
"""
function TestBase()
    rays = [-1 -1 -1; -1 -1 0; -1 -1 1; -1 -1 2; -1 -1 3; -1 -1 4; -1 -1 5; -1 0 -1; -1 0 0; -1 0 1; -1 0 2; -1 0 3; -1 0 4; -1 1 -1; -1 1 0; -1 1 1; -1 1 2; -1 1 3; -1 2 -1; -1 2 0; -1 2 1; -1 2 2; -1 3 -1; -1 3 0; -1 3 1; -1 4 -1; -1 4 0; -1 5 -1; 0 -1 -1; 0 -1 0; 0 -1 1; 0 -1 2; 0 0 -1; 0 0 1; 0 1 -1; 0 1 0; 0 2 -1; 1 -1 -1]
    cones = IncidenceMatrix([[36, 37, 38], [35, 37, 38], [34, 36, 38], [33, 35, 38], [32, 34, 38], [31, 32, 38], [30, 31, 38], [29, 33, 38], [29, 30, 38], [27, 28, 37], [26, 28, 37], [26, 27, 28], [25, 36, 37], [25, 27, 37], [24, 26, 27], [24, 25, 27], [23, 26, 37], [23, 24, 26], [22, 34, 36], [22, 25, 36], [21, 24, 25], [21, 22, 25], [20, 23, 24], [20, 21, 24], [19, 35, 37], [19, 23, 37], [19, 20, 23], [18, 32, 34], [18, 22, 34], [17, 21, 22], [17, 18, 22], [16, 20, 21], [16, 17, 21], [15, 19, 20], [15, 16, 20], [14, 33, 35], [14, 19, 35], [14, 15, 19], [13, 18, 32], [12, 17, 18], [12, 13, 18], [11, 16, 17], [11, 12, 17], [10, 15, 16], [10, 11, 16], [9, 14, 15], [9, 10, 15], [8, 29, 33], [8, 14, 33], [8, 9, 14], [7, 13, 32], [6, 12, 13], [6, 7, 32], [6, 7, 13], [5, 11, 12], [5, 6, 32], [5, 6, 12], [4, 31, 32], [4, 10, 11], [4, 5, 32], [4, 5, 11], [3, 30, 31], [3, 9, 10], [3, 4, 31], [3, 4, 10], [2, 29, 30], [2, 8, 9], [2, 3, 30], [2, 3, 9], [1, 8, 29], [1, 2, 29], [1, 2, 8]])
    return NormalToricVariety(PolyhedralFan(rays, cones))
end
export TestBase
