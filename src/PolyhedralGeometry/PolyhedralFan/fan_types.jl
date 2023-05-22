###############################################################################
###############################################################################
### Types for fan-like objects
###############################################################################
###############################################################################

# We introduce this abstract (hidden) type to allow for other objects to be
# used like polyhedral fans without duplicating too much code, concretely we
# want to be able to directly access rays, maximal_cones, etc for
# NormalToricVariety's.
# abstract type _FanLikeType{T} end
# abstract type AbstractNormalToricVariety end

struct PolyhedralFan{T}
    pm_fan::Polymake.BigObject
    parent_field::Field
    PolyhedralFan{T}(pm::Polymake.BigObject, f::Field) where T<:scalar_types = new{T}(pm, f)
    PolyhedralFan{QQFieldElem}(pm::Polymake.BigObject) = new{QQFieldElem}(pm, QQ)
end

const _FanLikeType = Union{NormalToricVarietyType, PolyhedralFan}
const _FanLikeTypeQQ = Union{NormalToricVarietyType, PolyhedralFan{QQFieldElem}}

get_scalar_type(::PolyhedralFan{T}) where T<:scalar_types = T
get_scalar_type(::NormalToricVarietyType) = QQFieldElem
coefficient_field(x::PolyhedralFan{T}) where T<:scalar_types = x.parent_field



################################################################################

# Lineality often causes certain collections to be empty;
# the following definition allows to easily construct a working empty SOI

_empty_access() = nothing

function _empty_subobjectiterator(::Type{T}, Obj:: _FanLikeType) where T
    return SubObjectIterator{T}(Obj, _empty_access, 0, NamedTuple())
end

for f in ("_point_matrix", "_vector_matrix", "_generator_matrix")
    M = Symbol(f)
    @eval begin
        function $M(::Val{_empty_access}, P::_FanLikeType; homogenized=false)
            scalar_regexp = match(r"[^<]*<(.*)>[^>]*", String(Polymake.type_name(pm_object(P))))
            typename = scalar_regexp[1]
            T = _scalar_type_to_polymake(scalar_type_to_oscar[typename])
            return Polymake.Matrix{T}(undef, 0, Polymake.polytope.ambient_dim(pm_object(P)) + homogenized)
        end
    end
end

for f in ("_facet_indices", "_ray_indices", "_vertex_indices", "_vertex_and_ray_indices")
    M = Symbol(f)
    @eval begin
        $M(::Val{_empty_access}, P::_FanLikeType) = return Polymake.IncidenceMatrix(0, Polymake.polytope.ambient_dim(P))
    end
end

for f in ("_linear_inequality_matrix", "_linear_equation_matrix")
    M = Symbol(f)
    @eval begin
        function $M(::Val{_empty_access}, P::_FanLikeType)
            scalar_regexp = match(r"[^<]*<(.*)>[^>]*", String(Polymake.type_name(P)))
            typename = scalar_regexp[1]
            T = _scalar_type_to_polymake(scalar_type_to_oscar[typename])
            return Polymake.Matrix{T}(undef, 0, Polymake.polytope.ambient_dim(P))
        end
    end
end

for f in ("_affine_inequality_matrix", "_affine_equation_matrix")
    M = Symbol(f)
    @eval begin
        function $M(::Val{_empty_access}, P::_FanLikeType)
            scalar_regexp = match(r"[^<]*<(.*)>[^>]*", String(Polymake.type_name(P)))
            typename = scalar_regexp[1]
            T = _scalar_type_to_polymake(scalar_type_to_oscar[typename])
            return Polymake.Matrix{T}(undef, 0, Polymake.polytope.ambient_dim(P) + 1)
        end
    end
end

_matrix_for_polymake(::Val{_empty_access}) = _point_matrix
