###
# Tropical varieties in Oscar
# ===========================
###



###
# 0. Definition
# -------------
# M = typeof(min) or typeof(max):
#   min or max convention, affecting initial ideals
# EMB = true or false:
#   embedded or abstract tropical variety
#   embedded tropical variety = weighted polyhedral complex in euclidean space
#   abstract tropical variety = weighted hypergraph with enumerated vertices
###

@attributes mutable struct TropicalVariety{M,EMB} <: TropicalVarietySupertype{M,EMB}
    polyhedralComplex::PolyhedralComplex
    function TropicalVariety{M,EMB}(Sigma::PolyhedralComplex) where {M,EMB}
        return new{M,EMB}(Sigma)
    end
end
export TropicalVariety

function pm_object(T::TropicalVariety)
    if has_attribute(T,:polymake_bigobject)
        return get_attribute(T,:polymake_bigobject)
    end
    error("pm_object(T::TropicalVariety): Has no polymake bigobject.")
end



###
# 1. Printing
# -----------
###

function Base.show(io::IO, tv::TropicalVariety{M, EMB}) where {M, EMB}
    if EMB
        print(io, "A $(repr(M)) tropical variety of dimension $(dim(tv)) embedded in $(ambient_dim(tv))-dimensional Euclidian space")
    else
        print(io, "An abstract $(repr(M)) tropical variety of dimension $(dim(tv))")
    end
end



###
# 2. Basic constructors
# ---------------------
###

@doc Markdown.doc"""
    TropicalVariety()

Construct the embedded tropical variety of a polynomial ideal over a (possibly trivially) valued field

# Examples
"""
# todo: Dartmouth
# function TropicalVariety()
#
#     return #...
# end



@doc Markdown.doc"""
    TropicalVariety{M,EMB}(Sigma::PolyhedralComplex)

Construct the abstract tropical variety from a polyhedral complex

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4]]);

julia> VR = [0 0; 1 0; 0 1; -1 -1];

julia> far_vertices = [2,3,4];

julia> Sigma = PolyhedralComplex(IM, VR, far_vertices);

julia> tropicalLine = TropicalVariety{min,true}(Sigma)
"""


###
# 3. Basic properties
# -------------------
###
