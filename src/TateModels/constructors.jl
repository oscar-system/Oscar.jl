##############################################
# 1: The Julia type for GlobalTateModel
##############################################

@attributes mutable struct GlobalTateModel
    a1::MPolyElem{fmpq}
    a2::MPolyElem{fmpq}
    a3::MPolyElem{fmpq}
    a4::MPolyElem{fmpq}
    a6::MPolyElem{fmpq}
    pt::MPolyElem{fmpq}
    base::Oscar.AbstractNormalToricVariety
    toric_ambient_space::Oscar.AbstractNormalToricVariety
    # TODO: Once new Oscar release is available
    #Y4::Oscar.ClosedSubvarietyOfToricVariety
    function GlobalTateModel(a1::MPolyElem{fmpq},
                             a2::MPolyElem{fmpq},
                             a3::MPolyElem{fmpq},
                             a4::MPolyElem{fmpq},
                             a6::MPolyElem{fmpq},
                             pt::MPolyElem{fmpq},
                             base::Oscar.AbstractNormalToricVariety,
                             toric_ambient_space::Oscar.AbstractNormalToricVariety)
                             # TODO: Once new Oscar release is available
                             #Y4::Oscar.ClosedSubvarietyOfToricVariety)
        return new(a1, a2, a3, a4, a6, pt, base, toric_ambient_space)
    end
end
export GlobalTateModel


#######################################
# 2: Generic constructor
#######################################

@doc Markdown.doc"""
    GenericGlobalTateModel(base::Oscar.AbstractNormalToricVariety)

This method constructs a global Tate model over a given toric base
3-fold. The Tate sections ``a_i`` are taken with (pseudo) random coefficients.

One way to achieve this is to first focus on the Cox ring of the toric
ambient space. This ring must be graded such that the Tate polynomial is
homogeneous and cuts out a Calabi-Yau hypersurface. Given this grading,
one can perform a triangulation. This triangulation will typically
take a long time to complete and yield a large number of candidate ambient
spaces. Typically, one wishes to focus on those spaces which contain the
base toric space in a manifest way. But even this criterion usually
allows for many ambient spaces.

This method here operates in the opposite way. It begins by extracting the rays and
maximal cones of the base space. Subsequently, those rays and cones are extended
to form one of the many toric 5-fold ambient spaces. The following example
demonstrates that this ambient space of a (singular) global Tate model need
not be smooth.

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

julia> t = GenericGlobalTateModel(base)
A global Tate model

julia> is_smooth(toric_ambient_space(t))
false
```
"""
function GenericGlobalTateModel(base::Oscar.AbstractNormalToricVariety)

    # Check if the base space is 3-dimensional
    if dim(base) != 3
        throw(ArgumentError("We currently focus on global Tate models over 3-dimensional toric base spaces"))
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

    # Compute the Tate sections
    f = hom(S1, S2, [gens(S2)[i] for i in 1:length(gens(S1))])
    a1 = sum([rand(Int) * f(b) for b in basis_of_global_sections(anticanonical_bundle(base))])
    a2 = sum([rand(Int) * f(b) for b in basis_of_global_sections(anticanonical_bundle(base)^2)])
    a3 = sum([rand(Int) * f(b) for b in basis_of_global_sections(anticanonical_bundle(base)^3)])
    a4 = sum([rand(Int) * f(b) for b in basis_of_global_sections(anticanonical_bundle(base)^4)])
    a6 = sum([rand(Int) * f(b) for b in basis_of_global_sections(anticanonical_bundle(base)^6)])

    # Compute the Tate polynomial
    x = gens(S2)[length(gens(S2))-2]
    y = gens(S2)[length(gens(S2))-1]
    z = gens(S2)[length(gens(S2))]
    pt = x^3 - y^2 + x*y*z*a1 + x^2*z^2*a2 + y*z^3*a3 + x*z^4*a4 + z^6*a6

    # TODO: Compute the toric hypersurface
    # TODO: Once new Oscar release is available
    #Y4 = Oscar.ClosedSubvarietyOfToricVariety(toric_ambient_space, [pt])

    # Return the global Tate model
    return GlobalTateModel(a1, a2, a3, a4, a6, pt, base, toric_ambient_space)

end
export GenericGlobalTateModel


@doc Markdown.doc"""
    GenericGlobalTateModelOverProjectiveSpace()

This method constructs a global Tate model over the 3-dimensional projective space.

# Examples
```jldoctest
julia> using Oscar

julia> GenericGlobalTateModelOverProjectiveSpace()
A global Tate model
```
"""
GenericGlobalTateModelOverProjectiveSpace() = GenericGlobalTateModel(projective_space(NormalToricVariety,3))
export GenericGlobalTateModelOverProjectiveSpace


@doc Markdown.doc"""
    SpecificGlobalTateModel(ais::Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}, base::Oscar.AbstractNormalToricVariety)

This method operates analogously to `GenericGlobalTateModel(base::Oscar.AbstractNormalToricVariety)`.
The only difference is that the Tate sections ``a_i`` can be specified with non-generic values.

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

julia> a1 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base))]);

julia> a2 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^2)]);

julia> a3 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^3)]);

julia> a4 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^4)]);

julia> a6 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)]);

julia> t = SpecificGlobalTateModel([a1, a2, a3, a4, a6], base)
A global Tate model

julia> is_smooth(toric_ambient_space(t))
false
```
"""
function SpecificGlobalTateModel(ais::Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}, base::Oscar.AbstractNormalToricVariety)

    # Extract information about the toric base
    base_rays = matrix(ZZ,rays(base))
    base_cones = matrix(ZZ, ray_indices(maximal_cones(base)))
    S1 = cox_ring(base)

    # Consistency checks
    if dim(base) != 3
        throw(ArgumentError("We currently focus on global Tate models over 3-dimensional toric base spaces"))
    end
    if length(ais) != 5
        throw(ArgumentError("We require exactly 5 Tate section"))
    end
    if any(k -> parent(k) != S1, ais)
        throw(ArgumentError("All Tate sections must reside in the Cox ring of the base toric variety"))
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

    # Compute the Tate polynomial
    x = gens(S2)[length(gens(S2))-2]
    y = gens(S2)[length(gens(S2))-1]
    z = gens(S2)[length(gens(S2))]
    ring_map = hom(S1, S2, [gens(S2)[i] for i in 1:length(gens(S1))])
    pt = x^3 - y^2 + x*y*z*ring_map(ais[1]) + x^2*z^2*ring_map(ais[2]) + y*z^3*ring_map(ais[3]) + x*z^4*ring_map(ais[4]) + z^6*ring_map(ais[5])

    # TODO: Compute the toric hypersurface
    # TODO: Once new Oscar release is available
    #Y4 = Oscar.ClosedSubvarietyOfToricVariety(toric_ambient_space, [pt])

    # Return the global Tate model
    return GlobalTateModel(ais[1], ais[2], ais[3], ais[4], ais[5], pt, base, toric_ambient_space)

end
export SpecificGlobalTateModel


#######################################
# 3: Display
#######################################

function Base.show(io::IO, cy::GlobalTateModel)
    join(io, "A global Tate model")
end
