######################################################
# 1: The Julia type for TateModelOverGeneralBaseSpace
######################################################

@attributes mutable struct TateModelOverGeneralBaseSpace
    a1::MPolyElem_dec{fmpq, fmpq_mpoly}
    a2::MPolyElem_dec{fmpq, fmpq_mpoly}
    a3::MPolyElem_dec{fmpq, fmpq_mpoly}
    a4::MPolyElem_dec{fmpq, fmpq_mpoly}
    a6::MPolyElem_dec{fmpq, fmpq_mpoly}
    pt::MPolyElem_dec{fmpq, fmpq_mpoly}
    auxiliary_base_space::Oscar.AbstractNormalToricVariety
    auxiliary_ambient_space::Oscar.AbstractNormalToricVariety
    # TODO: Once new Oscar release is available
    #Y4::Oscar.ClosedSubvarietyOfToricVariety
    function TateModelOverGeneralBaseSpace(a1::MPolyElem_dec{fmpq, fmpq_mpoly},
                                            a2::MPolyElem_dec{fmpq, fmpq_mpoly},
                                            a3::MPolyElem_dec{fmpq, fmpq_mpoly},
                                            a4::MPolyElem_dec{fmpq, fmpq_mpoly},
                                            a6::MPolyElem_dec{fmpq, fmpq_mpoly},
                                            pt::MPolyElem_dec{fmpq, fmpq_mpoly},
                                            auxiliary_base_space::Oscar.AbstractNormalToricVariety,
                                            auxiliary_ambient_space::Oscar.AbstractNormalToricVariety)
                                            # TODO: Once new Oscar release is available
                                            #Y4::Oscar.ClosedSubvarietyOfToricVariety)
        return new(a1, a2, a3, a4, a6, pt, auxiliary_base_space, auxiliary_ambient_space)
    end
end
export TateModelOverGeneralBaseSpace


#######################################
# 2: Generic constructor
#######################################

@doc Markdown.doc"""
    GlobalTateModel(ais::Vector{fmpq_mpoly}, auxiliary_base_ring::MPolyRing)

This method constructs a global Tate model over a base space that is not
fully specified. The following example exemplifies this approach.

# Examples
```jldoctest
julia> using Oscar

julia> auxiliary_base_ring, (a10,a21,a32,a43,a65,w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];

julia> a1 = a10;

julia> a2 = a21 * w;

julia> a3 = a32 * w^2;

julia> a4 = a43 * w^3;

julia> a6 = a65 * w^5;

julia> ais = [a1, a2, a3, a4, a6];

julia> t = GlobalTateModel(ais, auxiliary_base_ring)
A global Tate model over a not fully specified base

julia> tate_polynomial(t)
a10*x*y*z + a21*w*x^2*z^2 + a32*w^2*y*z^3 + a43*w^3*x*z^4 + a65*w^5*z^6 + x^3 - y^2

julia> auxiliary_base_space(t)
A normal, affine, 6-dimensional toric variety

julia> auxiliary_ambient_space(t)
A normal toric variety

julia> dim(auxiliary_ambient_space(t))
8
```
"""
function GlobalTateModel(ais::Vector{fmpq_mpoly}, auxiliary_base_ring::MPolyRing)

    # Check if the base space is 3-dimensional
    if length(ais) != 5
        throw(ArgumentError("We expect exactly 5 Tate sections"))
    end
    if any(k -> parent(k) != auxiliary_base_ring, ais)
        throw(ArgumentError("All Tate sections must reside in the provided auxiliary base ring"))
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
    f = hom(auxiliary_base_ring, S2, [gens(S2)[i] for i in 1:length(gens(auxiliary_base_ring))])
    x = gens(S2)[length(gens(S2))-2]
    y = gens(S2)[length(gens(S2))-1]
    z = gens(S2)[length(gens(S2))]
    pt = x^3 - y^2 + x*y*z*f(ais[1]) + x^2*z^2*f(ais[2]) + y*z^3*f(ais[3]) + x*z^4*f(ais[4]) + z^6*f(ais[5])

    # TODO: Compute the toric hypersurface
    # TODO: Once new Oscar release is available
    #Y4 = Oscar.ClosedSubvarietyOfToricVariety(auxiliary_ambient_space, [pt])

    # Return the global Tate model
    return TateModelOverGeneralBaseSpace(f(ais[1]), f(ais[2]), f(ais[3]), f(ais[4]), f(ais[5]), pt, auxiliary_base_space, auxiliary_ambient_space)

end
export GlobalTateModel


#######################################
# 3: Display
#######################################

function Base.show(io::IO, cy::TateModelOverGeneralBaseSpace)
    join(io, "A global Tate model over a not fully specified base")
end
