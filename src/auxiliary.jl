################################################################
# 1: Construct ambient space from given base
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
# 2: Construct the Weierstrass polynomial
################################################################

function _weierstrass_polynomial(base::Oscar.AbstractNormalToricVariety, toric_ambient_space::Oscar.AbstractNormalToricVariety)
    f = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^4)])
    g = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)])
    return _weierstrass_polynomial(base, toric_ambient_space, f, g)
end

function _weierstrass_polynomial(base::Oscar.AbstractNormalToricVariety, toric_ambient_space::Oscar.AbstractNormalToricVariety, f::MPolyElem_dec{fmpq, fmpq_mpoly}, g::MPolyElem_dec{fmpq, fmpq_mpoly})
    S = cox_ring(toric_ambient_space)
    x = gens(S)[length(gens(S))-2]
    y = gens(S)[length(gens(S))-1]
    z = gens(S)[length(gens(S))]
    ring_map = hom(parent(f), S, [gens(S)[i] for i in 1:length(gens(S))-3])
    return [ring_map(f), ring_map(g), x^3 - y^2 + ring_map(f)*x*z^4 + ring_map(g)*z^6]
end

function _weierstrass_polynomial(base::Oscar.AbstractNormalToricVariety, toric_ambient_space::Oscar.AbstractNormalToricVariety, f::fmpq_mpoly, g::fmpq_mpoly)
    S = cox_ring(toric_ambient_space)
    x = gens(S)[length(gens(S))-2]
    y = gens(S)[length(gens(S))-1]
    z = gens(S)[length(gens(S))]
    ring_map = hom(parent(f), S, [gens(S)[i] for i in 1:length(gens(S))-3])
    return [ring_map(f), ring_map(g), x^3 - y^2 + ring_map(f)*x*z^4 + ring_map(g)*z^6]
end


################################################################
# 3: Construct the Tate polynomial
################################################################

function _tate_polynomial(base::Oscar.AbstractNormalToricVariety, toric_ambient_space::Oscar.AbstractNormalToricVariety)
    a1 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base))])
    a2 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^2)])
    a3 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^3)])
    a4 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^4)])
    a6 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)])
    return _tate_polynomial(base, toric_ambient_space, [a1, a2, a3, a4, a6])
end

function _tate_polynomial(base::Oscar.AbstractNormalToricVariety, toric_ambient_space::Oscar.AbstractNormalToricVariety, ais::Vector{MPolyElem_dec{fmpq, fmpq_mpoly}})
    S = cox_ring(toric_ambient_space)
    x = gens(S)[length(gens(S))-2]
    y = gens(S)[length(gens(S))-1]
    z = gens(S)[length(gens(S))]
    ring_map = hom(parent(ais[1]), S, [gens(S)[i] for i in 1:length(gens(S))-3])
    (a1, a2, a3, a4, a6) = [ring_map(k) for k in ais]
    return [a1, a2, a3, a4, a6, x^3 - y^2 - x*y*z*a1 + x^2*z^2*a2 - y*z^3*a3 + x*z^4*a4 + z^6*a6]
end

function _tate_polynomial(base::Oscar.AbstractNormalToricVariety, toric_ambient_space::Oscar.AbstractNormalToricVariety, ais::Vector{fmpq_mpoly})
    S = cox_ring(toric_ambient_space)
    x = gens(S)[length(gens(S))-2]
    y = gens(S)[length(gens(S))-1]
    z = gens(S)[length(gens(S))]
    ring_map = hom(parent(ais[1]), S, [gens(S)[i] for i in 1:length(gens(S))-3])
    (a1, a2, a3, a4, a6) = [ring_map(k) for k in ais]
    return [a1, a2, a3, a4, a6, x^3 - y^2 - x*y*z*a1 + x^2*z^2*a2 - y*z^3*a3 + x*z^4*a4 + z^6*a6]
end
