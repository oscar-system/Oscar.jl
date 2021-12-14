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
    # GapTV::GapObj
    polymakeTV::Polymake.BigObject
end
export TropicalVariety



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
function TropicalVariety()
    # todo: Dartmouth
    return #...
end


@doc Markdown.doc"""
    TropicalVariety{M,EMB}()

Construct the abstract tropical variety from a list of vertices and maximal cells

# Examples
"""
function TropicalVariety{M, EMB}(pf::PolyhedralFan) where {M, EMB}
    if EMB
        return TropicalVariety{M, EMB}(polyhedral_complex_workaround(pm_object(pf)))
    else
        return TropicalVariety{M, EMB}(Polymake.fan.PolyhedralComplex(pm_object(pf)))
    end
end



###
# 3. Basic properties
# -------------------
###
