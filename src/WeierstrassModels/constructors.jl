##############################################
# 1: The Julia type for GlobalWeierstrassModel
##############################################

@attributes mutable struct GlobalWeierstrassModel
    poly_f::MPolyElem{fmpq}
    poly_g::MPolyElem{fmpq}
    pw::MPolyElem{fmpq}
    toric_base_space::Oscar.AbstractNormalToricVariety
    toric_ambient_space::Oscar.AbstractNormalToricVariety
    # TODO: Once new Oscar release is available
    #Y4::Oscar.ClosedSubvarietyOfToricVariety
    function GlobalWeierstrassModel(poly_f::MPolyElem{fmpq},
                             poly_g::MPolyElem{fmpq},
                             pw::MPolyElem{fmpq},
                             toric_base_space::Oscar.AbstractNormalToricVariety,
                             toric_ambient_space::Oscar.AbstractNormalToricVariety)
                             # TODO: Once new Oscar release is available
                             #Y4::Oscar.ClosedSubvarietyOfToricVariety)
        return new(poly_f, poly_g, pw, toric_base_space, toric_ambient_space)
    end
end
export GlobalWeierstrassModel


#######################################
# 2: Generic constructor
#######################################

@doc Markdown.doc"""
    GenericGlobalWeierstrassModel(base::Oscar.AbstractNormalToricVariety)

This method constructs a global Weierstrass model over a given toric base
3-fold. The Weierstrass sections ``f`` and ``g`` are taken with (pseudo)
random coefficients.

One way to achieve this is to first focus on the Cox ring of the toric
ambient space. This ring must be graded such that the Weierstrass polynomial
is homogeneous and cuts out a Calabi-Yau hypersurface. Given this grading,
one can perform a triangulation. This triangulation will typically
take a long time to complete and yield a large number of candidate ambient
spaces. Typically, one wishes to focus on those spaces which contain the
toric base space in a manifest way. But even this criterion usually
allows for many ambient spaces to be used.

To circumvent this obstacle, this method operates in the opposite direction, i.e.
"bottom-up". That is, it begins by extracting the rays and maximal cones of the chosen
toric base space. Subsequently, those rays and cones are extended
to form one of the many toric ambient spaces. The following example
demonstrates that the so-obtained ambient space of a (singular) global Weierstrass
model need not be smooth.

# Examples
```jldoctest
julia> using Oscar

julia> test_space = hirzebruch_surface(2) * projective_space(NormalToricVariety,1)
A normal toric variety

julia> test_space1 = blowup_on_ith_minimal_torus_orbit(test_space,1,"e1")
A normal toric variety

julia> test_space2 = blowup_on_ith_minimal_torus_orbit(test_space1,1,"e2")
A normal toric variety

julia> base = blowup_on_ith_minimal_torus_orbit(test_space2,1,"e3")
A normal toric variety

julia> w = GenericGlobalWeierstrassModel(base)
A global Weierstrass model over a concrete base

julia> is_smooth(toric_ambient_space(w))
false
```
"""
function GenericGlobalWeierstrassModel(base::Oscar.AbstractNormalToricVariety)

    # Check if the base space is 3-dimensional
    if dim(base) != 3
        throw(ArgumentError("We currently focus on global Weierstrass models over 3-dimensional toric base spaces"))
    end

    # Extract information about the toric base
    base_rays = matrix(ZZ,rays(base))
    base_cones = matrix(ZZ, ray_indices(maximal_cones(base)))
    S1 = cox_ring(base)

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

    # Construct the ambient space and perform consistency check
    toric_ambient_space = NormalToricVariety(PolyhedralFan(ambient_space_rays, ambient_space_max_cones; non_redundant = true))
    set_coordinate_names(toric_ambient_space, vcat(["u$i" for i in 1:nrays(base)], ["x", "y", "z"]))
    S2 = cox_ring(toric_ambient_space)

    # Compute the Weierstrass sections
    ring_map = hom(S1, S2, [gens(S2)[i] for i in 1:length(gens(S1))])
    f = sum([rand(Int) * ring_map(b) for b in basis_of_global_sections(anticanonical_bundle(base)^4)])
    g = sum([rand(Int) * ring_map(b) for b in basis_of_global_sections(anticanonical_bundle(base)^6)])

    # Compute the Weierstrass polynomial
    x = gens(S2)[length(gens(S2))-2]
    y = gens(S2)[length(gens(S2))-1]
    z = gens(S2)[length(gens(S2))]
    pw = x^3 - y^2 + f*x*z^4 + g*z^6

    # TODO: Compute the toric hypersurface
    # TODO: Once new Oscar release is available
    #Y4 = Oscar.ClosedSubvarietyOfToricVariety(toric_ambient_space, [pt])

    # Return the global Weierstrass model
    return GlobalWeierstrassModel(f, g, pw, base, toric_ambient_space)

end
export GenericGlobalWeierstrassModel


@doc Markdown.doc"""
    GenericGlobalWeierstrassModelOverProjectiveSpace()

This method constructs a global Weierstrass model over the 3-dimensional projective space.

# Examples
```jldoctest
julia> using Oscar

julia> GenericGlobalWeierstrassModelOverProjectiveSpace()
A global Weierstrass model over a concrete base
```
"""
GenericGlobalWeierstrassModelOverProjectiveSpace() = GenericGlobalWeierstrassModel(projective_space(NormalToricVariety,3))
export GenericGlobalWeierstrassModelOverProjectiveSpace



@doc Markdown.doc"""
    SpecificGlobalWeierstrassModel(f::MPolyElem_dec{fmpq, fmpq_mpoly}, g::MPolyElem_dec{fmpq, fmpq_mpoly}, base::Oscar.AbstractNormalToricVariety)

This method operates analogously to `GenericGlobalWeierstrassModel(base::Oscar.AbstractNormalToricVariety)`.
The only difference is that the Weierstrass sections ``f`` and ``g`` can be specified with non-generic values.

# Examples
```jldoctest
julia> using Oscar

julia> test_space = hirzebruch_surface(2) * projective_space(NormalToricVariety,1)
A normal toric variety

julia> test_space1 = blowup_on_ith_minimal_torus_orbit(test_space,1,"e1")
A normal toric variety

julia> test_space2 = blowup_on_ith_minimal_torus_orbit(test_space1,1,"e2")
A normal toric variety

julia> base = blowup_on_ith_minimal_torus_orbit(test_space2,1,"e3")
A normal toric variety

julia> f = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^4)]);

julia> g = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)]);

julia> w = SpecificGlobalWeierstrassModel(f, g, base)
A global Weierstrass model over a concrete base

julia> is_smooth(toric_ambient_space(w))
false
```
"""
function SpecificGlobalWeierstrassModel(f::MPolyElem_dec{fmpq, fmpq_mpoly}, g::MPolyElem_dec{fmpq, fmpq_mpoly}, base::Oscar.AbstractNormalToricVariety)

    # Extract information about the toric base
    base_rays = matrix(ZZ,rays(base))
    base_cones = matrix(ZZ, ray_indices(maximal_cones(base)))
    S1 = cox_ring(base)

    # Consistency checks
    if dim(base) != 3
        throw(ArgumentError("We currently focus on global Weierstrass models over 3-dimensional toric base spaces"))
    end
    if parent(f) != S1
        throw(ArgumentError("All Weierstrass sections must reside in the Cox ring of the base toric variety"))
    end
    if parent(g) != S1
        throw(ArgumentError("All Weierstrass sections must reside in the Cox ring of the base toric variety"))
    end

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

    # Construct the ambient space and perform consistency check
    toric_ambient_space = NormalToricVariety(PolyhedralFan(ambient_space_rays, ambient_space_max_cones; non_redundant = true))
    set_coordinate_names(toric_ambient_space, vcat(["u$i" for i in 1:nrays(base)], ["x", "y", "z"]))
    S2 = cox_ring(toric_ambient_space)

    # Compute the Weierstrass sections
    ring_map = hom(S1, S2, [gens(S2)[i] for i in 1:length(gens(S1))])

    # Compute the Weierstrass polynomial
    x = gens(S2)[length(gens(S2))-2]
    y = gens(S2)[length(gens(S2))-1]
    z = gens(S2)[length(gens(S2))]
    pw = x^3 - y^2 + ring_map(f)*x*z^4 + ring_map(g)*z^6

    # TODO: Compute the toric hypersurface
    # TODO: Once new Oscar release is available
    #Y4 = Oscar.ClosedSubvarietyOfToricVariety(toric_ambient_space, [pt])

    # Return the global Weierstrass model
    return GlobalWeierstrassModel(f, g, pw, base, toric_ambient_space)

end
export SpecificGlobalWeierstrassModel


#######################################
# 3: Display
#######################################

function Base.show(io::IO, cy::GlobalWeierstrassModel)
    join(io, "A global Weierstrass model over a concrete base")
end
