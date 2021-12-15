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
    function TropicalHypersurface{M,EMB}(pm::Polymake.BigObject) where {M,EMB}
        return new{M,EMB}(pm)
    end
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

Returns the tropical hypersurface of a tropical polynomial.

# Examples
```jldoctest
julia> T = tropical_numbers(min)
Tropical ring (min)

julia> Txy,(x,y) = T["x","y"]
(Multivariate Polynomial Ring in x, y over Tropical ring (min), AbstractAlgebra.Generic.MPoly{Oscar.TropicalRingElem{typeof(min)}}[x, y])

julia> f = x+y+1
x + y + (1)

julia> Tf = TropicalHypersurface(f)
A min tropical hypersurface embedded in 2-dimensional Euclidian space
```
"""
function TropicalHypersurface(f::Union{AbstractAlgebra.Generic.MPoly{Oscar.TropicalNumbersElem{typeof(min)}},
                                       AbstractAlgebra.Generic.MPoly{Oscar.TropicalNumbersElem{typeof(max)}}})
    if total_degree(f) <= 0
        error("Tropical hypersurfaces of constant polynomials not supported.")
    end
    convention = fun(base_ring(f))
    fstr = Tuple(tropical_polynomial_to_polymake(f))
    pmpoly = Polymake.common.totropicalpolynomial(fstr...)
    pmhyp = Polymake.tropical.Hypersurface{convention}(POLYNOMIAL=pmpoly)
    Vf = TropicalHypersurface{convention, true}(pmhyp)
    set_attribute!(Vf,:tropical_polynomial,f)
    return Vf
end

@doc Markdown.doc"""
    tropical_variety(f::Union{AbstractAlgebra.Generic.MPoly{Oscar.TropicalNumbersElem{typeof(min)}},
                              AbstractAlgebra.Generic.MPoly{Oscar.TropicalNumbersElem{typeof(max)}}})

Returns the tropical variety of a tropical polynomial in form of a TropicalHypersurface

# Examples
```jldoctest
julia> T = tropical_numbers(min)
Tropical ring (min)

julia> Txy,(x,y) = T["x","y"]
(Multivariate Polynomial Ring in x, y over Tropical ring (min), AbstractAlgebra.Generic.MPoly{Oscar.TropicalRingElem{typeof(min)}}[x, y])

julia> f = x+y+1
x + y + (1)

julia> Tf = TropicalHypersurface(f)
A min tropical hypersurface embedded in 2-dimensional Euclidian space
```
"""
function tropical_variety(f::Union{AbstractAlgebra.Generic.MPoly{Oscar.TropicalNumbersElem{typeof(min)}},
                                   AbstractAlgebra.Generic.MPoly{Oscar.TropicalNumbersElem{typeof(max)}}})
    return TropicalHypersurface(f)
end


@doc Markdown.doc"""
    TropicalHypersurface{M}(f::Union{AbstractAlgebra.Generic.MPoly{Oscar.TropicalNumbersElem{typeof(min)}},
                                     AbstractAlgebra.Generic.MPoly{Oscar.TropicalNumbersElem{typeof(max)}}})

Returns the tropical hypersurface of an algebraic polynomial.
If M=min, the tropical hypersurface will obey the min-convention.
If M=max, the tropical hypersurface will obey the max-convention.
If coefficient ring has a valuation, the tropical hypersurface will be constructed with respect to it.
If coefficient ring has no valuation, the tropical hypersurface will be constructed with respect to the trivial valuation.

# Examples
julia> K = PadicField(7, 2)

julia> Kxy, (x,y) = K["x", "y"]

julia> f = 7*x+y+49

julia> TropicalHypersurface{min}(f)

julia> TropicalHypersurface{max}(f)
"""
function TropicalHypersurface{M}(f::AbstractAlgebra.Generic.MPoly{<:RingElement}) where {M}
    tropf = tropical_polynomial(f,M)
    Tf = TropicalHypersurface(tropf)
    set_attribute!(Tf,:algebraic_polynomial,f)
    set_attribute!(Tf,:tropical_polynomial,tropf)
  return Tf
end


@doc Markdown.doc"""
    tropical_variety(f::AbstractAlgebra.Generic.MPoly{<:RingElement}, M::Union{typeof(min),typeof(max)})

Returns the tropical variety of an algebraic polynomial in the form of a TropicalHypersurface.
If M=min, the tropical hypersurface will obey the min-convention.
If M=max, the tropical hypersurface will obey the max-convention.
If coefficient ring has a valuation, the tropical hypersurface will be constructed with respect to it.
If coefficient ring has no valuation, the tropical hypersurface will be constructed with respect to the trivial valuation.
The function is the same as TropicalHypersurface{M}(f).

# Examples
julia> K = PadicField(7, 2)

julia> Kxy, (x,y) = K["x", "y"]

julia> f = 7*x+y+49

julia> tropical_variety(f,min)

julia> tropical_variety(f,max)
"""
function tropical_variety(f::AbstractAlgebra.Generic.MPoly{<:RingElement}, M::Union{typeof(min),typeof(max)})
    return TropicalHypersurface{M}(f)
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
