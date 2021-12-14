###
# Tropical curves in Oscar
# ========================
###



###
# 1. Definition
# -------------
# M = typeof(min) or typeof(max):
#   min or max convention, affecting initial ideals
# EMB = true or false:
#   embedded or abstract tropical curves
#   embedded tropical variety = graph embedded in euclidean space with weighted edges and vertices
#   abstract tropical variety = graph with enumerated vertices with weighted edges and vertices
###

@attributes mutable struct TropicalCurve{M,T} <: TropicalVarietySupertype{M,T}
    polymakeTV::Polymake.BigObject
    algebraicTV
end
export TropicalCurve
function pm_object(v::TropicalCurve)
  return v.polymakeTV
end



###
# 2. Basic constructors
# ---------------------
###

@doc Markdown.doc"""
    TropicalCurve{M, EMB}()

Construct an abstract of embedded tropical curve from a list of edges and a
vector of their lengths.

# Examples
"""
function TropicalCurve{M, EMB}() where {M, EMB}
    if EMB
        return #...
    else
        return #...
    end
end



###
# 3. Basic properties
# -------------------
###
