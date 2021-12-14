###
# Tropical hypersurfaces in Oscar
# ===============================
###



###
# 0. Definition
# -------------
# M = typeof(min) or typeof(max):
#   min or max convention, affecting initial forms, tropical polynomial, etc.
# EMB = true or false:
#   embedded or abstract tropical hypersurface
#   embedded tropical variety = weighted polyhedral complex in euclidean space
#   abstract tropical variety = weighted hypergraph with enumerated vertices
###

@attributes mutable struct TropicalHypersurface{M,EMB} <: TropicalVarietySupertype{M,EMB}
    # GapTV::GapObj
    polymakeTV::Polymake.BigObject
end
export TropicalHypersurface



###
# 1. Printing
# -----------
###

function Base.show(io::IO, th::TropicalHypersurface{M, EMB}) where {M, EMB}
    if EMB
        print(io, "A $(repr(M)) tropical hypersurface embedded in $(ambient_dim(th))-dimensional Euclidian space")
    else
        print(io, "An abstract $(repr(M)) tropical hypersurface of dimension $(dim(th))")
    end
end



###
# 2. Basic constructors
# ---------------------
###

@doc Markdown.doc"""
    TropicalHypersurface()

Construct the tropical hypersurface of a polynomial over a valued field

# Examples
"""
function TropicalHypersurface()
  return #...
end


@doc Markdown.doc"""
    TropicalHypersurface()

Construct the tropical hypersurface of a polynomial over the tropical numbers

# Examples
```jldoctest
julia> T = tropical_ring(min)
Tropical ring (min)

julia> Txy,(x,y) = T["x","y"]
(Multivariate Polynomial Ring in x, y over Tropical ring (min), AbstractAlgebra.Generic.MPoly{Oscar.TropicalRingElem{typeof(min)}}[x, y])

julia> f = x+y+1
x + y + (1)

julia> hyp = TropicalHypersurface(f)
A min tropical hypersurface embedded in 2-dimensional Euclidian space
```
"""
function TropicalHypersurface(f)
    if total_degree(f) <= 0
        error("Tropical variety of constant polynomials not supported.")
    end
    convention = fun(base_ring(f))
    fstr = Tuple(tropical_polynomial_to_polymake(f))
    pmpoly = Polymake.common.totropicalpolynomial(fstr...)
    pmhyp = Polymake.tropical.Hypersurface{convention}(POLYNOMIAL=pmpoly)
    return TropicalHypersurface{convention, true}(pmhyp)
end



###
# 3. Basic properties
# -------------------
###

@doc Markdown.doc"""
    dual_subdivision(TH::TropicalHypersurface{M, EMB})

Returns the dual subdivision of `TH` if it is embedded. Returns error otherwise

# Examples
A tropical hypersurface in RR^n is always of dimension n-1
```jldoctest
julia> T = tropical_numbers(min);

julia> Txy,(x,y) = T["x","y"];

julia> f = x+y+1;

julia> tropicalLine = TropicalHypersurface(f);

julia> dual_subdivision(tropicalLine)
# todo: add examples for varieties, curves and linear spaces
```
"""
function dual_subdivision(TH::TropicalHypersurface{M,EMB}) where {M,EMB}
    # not sure whether it makes sense to support abstract tropical hypersurfaces, but it can't hurt to check
    if !EMB
        error("tropical hypersurface not embedded")
    end

    return SubdivisionOfPoints(pm_object(TH).DUAL_SUBDIVISION)
end
export dual_subdivision
