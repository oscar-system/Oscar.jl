################################################################################
#
#  Abstract tropical variety supertype
#  ===================================
#
#  Not for public use, see concrete subtypes:
#  - TropicalVariety, see variety,jl
#  - TropicalHypersurface, see hypersurface.jl
#  - TropicalCurve, see curve.jl
#  - TropicalLinearSpace, see linear_space.jl
#
################################################################################

# minOrMax distinguishes between min- and max-convention
# isEmbedded distinguishes embedded and abstract tropical varieties
abstract type TropicalVarietySupertype{minOrMax,isEmbedded} end



################################################################################
#
#  Properties
#
################################################################################
# if TropV isa TropicalLinearSpace, then polymake_object isa Polymake.matroid.ValuatedMatroid
# if TropV isa TropicalHypersurface, then polymake_object isa Polymake.tropical.Hypersurface
function pm_object(TropV::TropicalVarietySupertype)
    @req has_attribute(TropV,:polymake_object) "no polymake object cached"
    return get_attribute(TropV,:polymake_object)
end


@doc raw"""
    polyhedral_complex(TropV::TropicalVariety)

Return the polyhedral complex of a tropical variety.
"""
function polyhedral_complex(TV::TropicalVarietySupertype)
    return TV.polyhedralComplex
end


@doc raw"""
    ambient_dim(TropV::TropicalVariety)

Return the ambient dimension of `TropV`.  Requires `TropV` to be embedded.
"""
function ambient_dim(TropV::TropicalVarietySupertype{minOrMax,true}) where minOrMax
    return ambient_dim(TropV.polyhedralComplex)
end


@doc raw"""
    codim(TropV::TropicalVariety)

Return the codimension of `TropV`.  Requires `TropV` to be embedded.
"""
function codim(TropV::TropicalVarietySupertype{minOrMax,true}) where minOrMax
    return codim(TropV.polyhedralComplex)
end


function convention(TropV::TropicalVarietySupertype{typeof(min),isEmbedded}) where isEmbedded
    return min
end
function convention(TropV::TropicalVarietySupertype{typeof(max),isEmbedded}) where isEmbedded
    return max
end


@doc raw"""
    dim(TropV::TropicalVariety)

Return the dimension of `TropV`.
"""
function dim(TropV::TropicalVarietySupertype)
    return dim(TropV.polyhedralComplex)
end


@doc raw"""
    f_vector(TropV::TropicalVariety)

Return the f-Vector of `TropV`.
"""
function f_vector(TropV::TropicalVarietySupertype)
    return f_vector(TropV.polyhedralComplex)
end


@doc raw"""
    lineality_dim(TropV::TropicalVariety)

Return the dimension of the lineality space of `TropV`.
"""
function lineality_dim(TropV::TropicalVarietySupertype)
    return lineality_dim(TropV.polyhedralComplex)
end


@doc raw"""
    lineality_space(TropV::TropicalVariety)

Return the lineality space of `TropV`.
"""
function lineality_space(TropV::TropicalVarietySupertype)
    return lineality_space(TropV.polyhedralComplex)
end


@doc raw"""
    maximal_polyhedra(TropV::TropicalVariety)

Return the maximal polyhedra of `TropV`.
"""
function maximal_polyhedra(TropV::TropicalVarietySupertype)
    return maximal_polyhedra(TropV.polyhedralComplex)
end


@doc raw"""
    maximal_polyhedra_and_multiplicities(TropV::TropicalVariety)

Return the maximal polyhedra of `TropV`.
"""
function maximal_polyhedra_and_multiplicities(TropV::TropicalVarietySupertype)
    TropVmults = multiplicities(TropV)
    return collect(zip(maximal_polyhedra(TropV),multiplicities(TropV)))
end


@doc raw"""
    minimal_faces(TropV::TropicalVariety)

Return the minimal faces of `TropV`.
"""
function minimal_faces(TropV::TropicalVarietySupertype)
    return minimal_faces(TropV.polyhedralComplex)
end


@doc raw"""
    multiplicities(TropV::TropicalVariety)

Return the multiplicities of `TropV`.
"""
function multiplicities(TropV::TropicalVarietySupertype)
    return TropV.multiplicities
end


@doc raw"""
    n_maximal_polyhedra(TropV::TropicalVariety)

Return the number of maximal polyhedra of `TropV`.
"""
function n_maximal_polyhedra(TropV::TropicalVarietySupertype)
    return n_maximal_polyhedra(TropV.polyhedralComplex)
end


@doc raw"""
    n_polyhedra(TropV::TropicalVariety)

Return the number of polyhedra of `TropV`.
"""
function n_polyhedra(TropV::TropicalVarietySupertype)
    return n_polyhedra(TropV.polyhedralComplex)
end


@doc raw"""
    n_vertices(TropV::TropicalVariety)

Return the number of vertices of `TropV`.
"""
function n_vertices(TropV::TropicalVarietySupertype)
    return n_vertices(TropV.polyhedralComplex)
end


@doc raw"""
    is_pure(TropV::TropicalVariety)

Return `true` if `TropV` is pure as a polyhedral complex, `false` otherwise.
"""
function is_pure(TropV::TropicalVarietySupertype)
    return is_pure(TropV.polyhedralComplex)
end


@doc raw"""
    is_simplicial(TropV::TropicalVariety)

Return `true` if `TropV` is simplicial as a polyhedral complex, `false` otherwise.
"""
function is_simplicial(TropV::TropicalVarietySupertype)
    return is_simplicial(TropV.polyhedralComplex)
end


@doc raw"""
    rays(TropV::TropicalVariety)

Return the rays of `TropV`.
"""
function rays(as::Type{RayVector{S}}, TropV::TropicalVarietySupertype{minOrMax,isEmbedded}) where {S,minOrMax,isEmbedded}
    return rays(as,TropV.polyhedralComplex)
end
function rays(TropV::TropicalVarietySupertype)
    return rays(TropV.polyhedralComplex)
end


@doc raw"""
    rays_modulo_lineality(TropV::TropicalVariety)

Return the rays modulo the lineality space of `TropV`.
"""
function rays_modulo_lineality(as::Type{RayVector{S}}, TropV::TropicalVarietySupertype{minOrMax,isEmbedded}) where {S,minOrMax,isEmbedded}
    return rays_modulo_lineality(as,TropV.polyhedralComplex)
end
function rays_modulo_lineality(TropV::TropicalVarietySupertype)
    return rays_modulo_lineality(TropV.polyhedralComplex)
end


@doc raw"""
    vertices_and_rays(TropV::TropicalVariety)

Return the vertices and rays of `TropV`.
"""
function vertices_and_rays(TropV::TropicalVarietySupertype)
    return vertices_and_rays(TropV.polyhedralComplex)
end


@doc raw"""
    vertices(TropV::TropicalVariety)

Return the vertices of `TropV`.
"""
function vertices(as::Type{PointVector{S}}, TropV::TropicalVarietySupertype{minOrMax,isEmbedded}) where {S,minOrMax,isEmbedded}
    return vertices(as,TropV.polyhedralComplex)
end
function vertices(TropV::TropicalVarietySupertype{minOrMax, isEmbedded}) where {minOrMax,isEmbedded}
    return vertices(TropV.polyhedralComplex)
end


@doc raw"""
    visualize(TropV::TropicalVariety)

Visualize `TropV`.
"""
function visualize(TropV::TropicalVarietySupertype)
    return visualize(TropV.polyhedralComplex)
end
