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

See [`ambient_dim(::PolyhedralComplex)`](@ref ambient_dim(PC::PolyhedralComplex)).
"""
function ambient_dim(TropV::TropicalVarietySupertype{minOrMax,true}) where minOrMax
    return ambient_dim(TropV.polyhedralComplex)
end


@doc raw"""
    codim(TropV::TropicalVariety)

See [`codim(::PolyhedralComplex)`](@ref codim(PC::PolyhedralComplex)).
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

See [`dim(::PolyhedralComplex)`](@ref dim(PC::PolyhedralComplex)).
"""
function dim(TropV::TropicalVarietySupertype)
    return dim(TropV.polyhedralComplex)
end


@doc raw"""
    f_vector(TropV::TropicalVariety)

See [`f_vector(::PolyhedralComplex)`](@ref f_vector(PC::PolyhedralComplex)).
"""
function f_vector(TropV::TropicalVarietySupertype)
    return f_vector(TropV.polyhedralComplex)
end


@doc raw"""
    lineality_dim(TropV::TropicalVariety)

See [`lineality_dim(::PolyhedralComplex)`](@ref lineality_dim(PC::PolyhedralComplex)).
"""
function lineality_dim(TropV::TropicalVarietySupertype)
    return lineality_dim(TropV.polyhedralComplex)
end


@doc raw"""
    lineality_space(TropV::TropicalVariety)

See [`lineality_space(::PolyhedralComplex)`](@ref lineality_space(PC::PolyhedralComplex{T}) where {T<:scalar_types}).
"""
function lineality_space(TropV::TropicalVarietySupertype)
    return lineality_space(TropV.polyhedralComplex)
end


@doc raw"""
    maximal_polyhedra(TropV::TropicalVariety)

See [`maximal_polyhedra(::PolyhedralComplex)`](@ref maximal_polyhedra(PC::PolyhedralComplex{T}) where {T<:scalar_types}).
"""
function maximal_polyhedra(TropV::TropicalVarietySupertype)
    return maximal_polyhedra(TropV.polyhedralComplex)
end


@doc raw"""
    maximal_polyhedra_and_multiplicities(TropV::TropicalVariety)

Return the maximal polyhedra and multiplicities of `TropV`.

# Examples
```jldoctest
julia> R,(x1,x2) = polynomial_ring(QQ,4);

julia> nu = tropical_semiring_map(QQ,2);

julia> f = 2*x1^2+x1*x2+x2^2+1
2*x1^2 + x1*x2 + x2^2 + 1

julia> TropH = tropical_hypersurface(f,nu)
Min tropical hypersurface

julia> maximal_polyhedra_and_multiplicities(TropH)
5-element Vector{Tuple{Polyhedron{QQFieldElem}, ZZRingElem}}:
 (Polyhedron in ambient dimension 4, 1)
 (Polyhedron in ambient dimension 4, 1)
 (Polyhedron in ambient dimension 4, 2)
 (Polyhedron in ambient dimension 4, 1)
 (Polyhedron in ambient dimension 4, 2)

```
"""
function maximal_polyhedra_and_multiplicities(TropV::TropicalVarietySupertype)
    TropVmults = multiplicities(TropV)
    return collect(zip(maximal_polyhedra(TropV),multiplicities(TropV)))
end


@doc raw"""
    minimal_faces(TropV::TropicalVariety)

See [`minimal_faces(::PolyhedralComplex)`](@ref minimal_faces(PC::PolyhedralComplex{T})  where {T<:scalar_types}).
"""
function minimal_faces(TropV::TropicalVarietySupertype)
    return minimal_faces(TropV.polyhedralComplex)
end


@doc raw"""
    multiplicities(TropV::TropicalVariety)

Return the multiplicities of `TropV`.  Order is the same as `maximal_polyhedra`.

# Examples
```jldoctest
julia> R,(x1,x2) = polynomial_ring(QQ,4);

julia> nu = tropical_semiring_map(QQ,2);

julia> f = 2*x1^2+x1*x2+x2^2+1
2*x1^2 + x1*x2 + x2^2 + 1

julia> TropH = tropical_hypersurface(f,nu)
Min tropical hypersurface

julia> multiplicities(TropH)
5-element Vector{ZZRingElem}:
 1
 1
 2
 1
 2

```
"""
function multiplicities(TropV::TropicalVarietySupertype)
    return TropV.multiplicities
end


@doc raw"""
    n_maximal_polyhedra(TropV::TropicalVariety)

See [`n_maximal_polyhedra(::PolyhedralComplex)`](@ref n_maximal_polyhedra(PC::PolyhedralComplex)).
"""
function n_maximal_polyhedra(TropV::TropicalVarietySupertype)
    return n_maximal_polyhedra(TropV.polyhedralComplex)
end


@doc raw"""
    n_polyhedra(TropV::TropicalVariety)

See [`n_polyhedra(::PolyhedralComplex)`](@ref n_polyhedra(PC::PolyhedralComplex)).
"""
function n_polyhedra(TropV::TropicalVarietySupertype)
    return n_polyhedra(TropV.polyhedralComplex)
end


@doc raw"""
    n_vertices(TropV::TropicalVariety)

See [`n_vertices(::PolyhedralComplex)`](@ref n_vertices(PC::PolyhedralComplex)).
"""
function n_vertices(TropV::TropicalVarietySupertype)
    return n_vertices(TropV.polyhedralComplex)
end


@doc raw"""
    is_pure(TropV::TropicalVariety)

See [`is_pure(::PolyhedralComplex)`](@ref is_pure(PC::PolyhedralComplex)).
"""
function is_pure(TropV::TropicalVarietySupertype)
    return is_pure(TropV.polyhedralComplex)
end


@doc raw"""
    is_simplicial(TropV::TropicalVariety)

See [`is_simplicial(::PolyhedralComplex)`](@ref is_simplicial(PC::PolyhedralComplex)).
"""
function is_simplicial(TropV::TropicalVarietySupertype)
    return is_simplicial(TropV.polyhedralComplex)
end


@doc raw"""
    rays(TropV::TropicalVariety)

See [`rays(::PolyhedralComplex)`](@ref rays(PC::PolyhedralComplex{T}) where {T<:scalar_types}).
"""
function rays(as::Type{RayVector{S}}, TropV::TropicalVarietySupertype{minOrMax,isEmbedded}) where {S,minOrMax,isEmbedded}
    return rays(as,TropV.polyhedralComplex)
end
function rays(TropV::TropicalVarietySupertype)
    return rays(TropV.polyhedralComplex)
end


@doc raw"""
    rays_modulo_lineality(TropV::TropicalVariety)

See [`rays_modulo_lineality(::PolyhedralComplex)`](@ref rays_modulo_lineality(PC::PolyhedralComplex{T}) where {T<:scalar_types}).
"""
function rays_modulo_lineality(as::Type{RayVector{S}}, TropV::TropicalVarietySupertype{minOrMax,isEmbedded}) where {S,minOrMax,isEmbedded}
    return rays_modulo_lineality(as,TropV.polyhedralComplex)
end
function rays_modulo_lineality(TropV::TropicalVarietySupertype)
    return rays_modulo_lineality(TropV.polyhedralComplex)
end


@doc raw"""
    vertices_and_rays(TropV::TropicalVariety)

See [`vertices_and_rays(::PolyhedralComplex)`](@ref vertices_and_rays(PC::PolyhedralComplex{T}) where {T<:scalar_types}).
"""
function vertices_and_rays(TropV::TropicalVarietySupertype)
    return vertices_and_rays(TropV.polyhedralComplex)
end


@doc raw"""
    vertices(TropV::TropicalVariety)

See [`vertices(::PolyhedralComplex)`](@ref vertices(as::Type{PointVector{T}}, PC::PolyhedralComplex{T}) where {T<:scalar_types}).
"""
function vertices(as::Type{PointVector{S}}, TropV::TropicalVarietySupertype{minOrMax,isEmbedded}) where {S,minOrMax,isEmbedded}
    return vertices(as,TropV.polyhedralComplex)
end
function vertices(TropV::TropicalVarietySupertype{minOrMax, isEmbedded}) where {minOrMax,isEmbedded}
    return vertices(TropV.polyhedralComplex)
end


@doc raw"""
    visualize(TropV::TropicalVariety)

See [`visualize(::PolyhedralComplex)`](@ref).
"""
function visualize(TropV::TropicalVarietySupertype)
    return visualize(TropV.polyhedralComplex)
end
