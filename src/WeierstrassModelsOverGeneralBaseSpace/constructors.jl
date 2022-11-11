######################################################
# 1: The Julia type for WeierstrassModelOverGeneralBaseSpace
######################################################

@attributes mutable struct WeierstrassModelOverGeneralBaseSpace
    f::MPolyElem_dec{fmpq, fmpq_mpoly}
    g::MPolyElem_dec{fmpq, fmpq_mpoly}
    pw::MPolyElem_dec{fmpq, fmpq_mpoly}
    auxiliary_base_space::Oscar.AbstractNormalToricVariety
    auxiliary_ambient_space::Oscar.AbstractNormalToricVariety
    # TODO: Once new Oscar release is available
    #Y4::Oscar.ClosedSubvarietyOfToricVariety
    function WeierstrassModelOverGeneralBaseSpace(f::MPolyElem_dec{fmpq, fmpq_mpoly},
                                                  g::MPolyElem_dec{fmpq, fmpq_mpoly},
                                                  pw::MPolyElem_dec{fmpq, fmpq_mpoly},
                                                  auxiliary_base_space::Oscar.AbstractNormalToricVariety,
                                                  auxiliary_ambient_space::Oscar.AbstractNormalToricVariety)
                                                  # TODO: Once new Oscar release is available
                                                  #Y4::Oscar.ClosedSubvarietyOfToricVariety)
        return new(f, g, pw, auxiliary_base_space, auxiliary_ambient_space)
    end
end
export WeierstrassModelOverGeneralBaseSpace


#######################################
# 2: Generic constructor
#######################################

@doc Markdown.doc"""
    GlobalWeierstrassModel(f::fmpq_mpoly, g::fmpq_mpoly, auxiliary_base_ring::MPolyRing)

This method constructs a global Weierstrass model over a base space that is not
fully specified. The following example illustrates this approach.

# Examples
```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (f, g) = QQ["f", "g"];

julia> w = GlobalWeierstrassModel(f, g, auxiliary_base_ring)
A global Weierstrass model over a general base space

julia> weierstrass_polynomial(w)
f*x*z^4 + g*z^6 + x^3 - y^2

julia> auxiliary_base_space(w)
A normal, affine, 2-dimensional toric variety

julia> auxiliary_ambient_space(w)
A normal toric variety

julia> dim(auxiliary_ambient_space(w))
4
```
"""
function GlobalWeierstrassModel(f::fmpq_mpoly, g::fmpq_mpoly, auxiliary_base_ring::MPolyRing)

    # Check if the base space is 3-dimensional
    if (parent(f) != auxiliary_base_ring) || (parent(g) != auxiliary_base_ring)
        throw(ArgumentError("All Weierstrass sections must reside in the provided auxiliary base ring"))
    end

    # Construct auxiliary base space
    auxiliary_base_space = affine_space(NormalToricVariety, length(gens(auxiliary_base_ring)))
    set_coordinate_names(auxiliary_base_space, [string(k) for k in gens(auxiliary_base_ring)])

    # Extract information about the auxiliary toric base space
    base_rays = matrix(ZZ,rays(auxiliary_base_space))
    base_cones = matrix(ZZ, ray_indices(maximal_cones(auxiliary_base_space)))

    # Construct the rays for the fibre ambient space
    xray = [0 for i in 1:ncols(base_rays)+2]
    yray = [0 for i in 1:ncols(base_rays)+2]
    zray = [0 for i in 1:ncols(base_rays)+2]
    xray[ncols(base_rays)+1] = 1
    yray[ncols(base_rays)+2] = 1
    zray[ncols(base_rays)+1] = -2
    zray[ncols(base_rays)+2] = -3

    # Construct the rays of the auxiliary toric ambient space
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
    auxiliary_ambient_space = NormalToricVariety(PolyhedralFan(ambient_space_rays, ambient_space_max_cones; non_redundant = true))
    set_coordinate_names(auxiliary_ambient_space, vcat([string(k) for k in gens(auxiliary_base_ring)], ["x", "y", "z"]))
    S2 = cox_ring(auxiliary_ambient_space)

    # Compute the Tate polynomial
    ring_map = hom(auxiliary_base_ring, S2, [gens(S2)[i] for i in 1:length(gens(auxiliary_base_ring))])
    x = gens(S2)[length(gens(S2))-2]
    y = gens(S2)[length(gens(S2))-1]
    z = gens(S2)[length(gens(S2))]
    pw = x^3 - y^2 + ring_map(f)*x*z^4 + ring_map(g)*z^6

    # TODO: Compute the toric hypersurface
    # TODO: Once new Oscar release is available
    #Y4 = Oscar.ClosedSubvarietyOfToricVariety(auxiliary_ambient_space, [pt])

    # Return the global Weierstrass model
    return WeierstrassModelOverGeneralBaseSpace(ring_map(f), ring_map(g), pw, auxiliary_base_space, auxiliary_ambient_space)

end
export GlobalWeierstrassModel


#######################################
# 3: Display
#######################################

function Base.show(io::IO, w::WeierstrassModelOverGeneralBaseSpace)
    join(io, "A global Weierstrass model over a general base space")
end
