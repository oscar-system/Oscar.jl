###
# Tropical hypersurfaces in Oscar
###

export TropicalHypersurface,
       dualsubdivision,
       polynomial

################################################################################
#
#  Definition
#
################################################################################

# We use M to record whether we are in the min/max case
# M is either typeof(min) or typeof(max)
# We use EMB to record whether the hypersurface is embedded or abstract
# EMB is either true or false:
#   embedded tropical variety = weighted polyhedral complex in euclidean space
#   abstract tropical variety = weighted hypergraph with enumerated vertices

@attributes mutable struct TropicalHypersurface{M,EMB} <: TropicalVarietySupertype{M,EMB}
    polyhedralComplex::PolyhedralComplex
    function TropicalHypersurface{M,EMB}(Sigma::PolyhedralComplex) where {M,EMB}
        if codim(Sigma)!=1
            error("TropicalHypersurface: input polyhedral complex not one-codimensional")
        end
        return new{M,EMB}(Sigma)
    end
end

function pm_object(T::TropicalHypersurface)
    if has_attribute(T,:polymake_bigobject)
        return get_attribute(T,:polymake_bigobject)
    end
    error("pm_object(T::TropicalHypersurface): Has no polymake bigobject")
end

################################################################################
#
#  Printing
#
################################################################################


function Base.show(io::IO, th::TropicalHypersurface{M, EMB}) where {M, EMB}
    if EMB
        print(io, "A $(repr(M)) tropical hypersurface embedded in $(ambient_dim(th))-dimensional Euclidian space")
    else
        print(io, "An abstract $(repr(M)) tropical hypersurface of dimension $(dim(th))")
    end
end

################################################################################
#
#  Constructors for tropical hypersurfaces
#
################################################################################

@doc Markdown.doc"""
    TropicalHypersurface(f::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{T}})

Return the tropical hypersurface defined by a tropical polynomial.

# Examples
```jldoctest
julia> T = TropicalSemiring(min)
Tropical semiring (min)

julia> Txy,(x,y) = T["x","y"]
(Multivariate Polynomial Ring in x, y over Tropical semiring (min), AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}}[x, y])

julia> f = x+y+1
x + y + (1)

julia> Tf = TropicalHypersurface(f)
A min tropical hypersurface embedded in 2-dimensional Euclidian space
```
"""
function TropicalHypersurface(f::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{T}}) where T
    if total_degree(f) <= 0
        error("Tropical hypersurfaces of constant polynomials not supported.")
    end
    M = convention(base_ring(f))
    fstr = Tuple(tropical_polynomial_to_polymake(f))
    pmpoly = Polymake.common.totropicalpolynomial(fstr...)
    pmhypproj = Polymake.tropical.Hypersurface{M}(POLYNOMIAL=pmpoly)
    pmhyp = Polymake.tropical.affine_chart(pmhypproj)
    Vf = TropicalHypersurface{M, true}(PolyhedralComplex(pmhyp))
    w = pmhypproj.WEIGHTS
    set_attribute!(Vf,:polymake_bigobject,pmhypproj)
    set_attribute!(Vf,:tropical_polynomial,f)
    set_attribute!(Vf,:weights,w)
    return Vf
end

# @doc Markdown.doc"""
#     tropical_variety(f::Union{AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}},
#                               AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(max)}}})

# Return the tropical variety defined by a tropical polynomial in form of a TropicalHypersurface

# # Examples
# ```jldoctest
# julia> T = TropicalSemiring(min)
# Tropical ring (min)

# julia> Txy,(x,y) = T["x","y"]
# (Multivariate Polynomial Ring in x, y over Tropical ring (min), AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}}[x, y])

# julia> f = x+y+1
# x + y + (1)

# julia> Tf = TropicalHypersurface(f)
# A min tropical hypersurface embedded in 2-dimensional Euclidian space
# ```
# """
# function tropical_variety(f::Union{AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}},
#                                    AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(max)}}})
#     return TropicalHypersurface(f)
# end

@doc Markdown.doc"""
    TropicalHypersurface{M}(f::MPolyElem,M::Union{typeof(min),typeof(max)}=min)

Return the tropical hypersurface defined by the tropicalization of an algebraic polynomial.
If M=min, the tropical hypersurface will obey the min-convention.
If M=max, the tropical hypersurface will obey the max-convention.
If coefficient ring has a valuation, the tropical hypersurface will be constructed with respect to it.
If coefficient ring has no valuation, the tropical hypersurface will be constructed with respect to the trivial valuation.

# Examples
```jldoctest
julia> K = PadicField(7, 2);

julia> Kxy, (x,y) = K["x", "y"]
(Multivariate Polynomial Ring in x, y over Field of 7-adic numbers, AbstractAlgebra.Generic.MPoly{padic}[x, y])

julia> f = 7*x+y+49;

julia> TropicalHypersurface(f, min)
A min tropical hypersurface embedded in 2-dimensional Euclidian space

julia> TropicalHypersurface(f, max)
A max tropical hypersurface embedded in 2-dimensional Euclidian space
```
"""
function TropicalHypersurface(f::MPolyElem,M::Union{typeof(min),typeof(max)}=min)
    tropf = tropical_polynomial(f,M)
    Tf = TropicalHypersurface(tropf)
    w = pm_object(Tf).WEIGHTS
    set_attribute!(Tf,:algebraic_polynomial,f)
    set_attribute!(Tf,:tropical_polynomial,tropf)
    set_attribute!(Tf,:weights,w)
    return Tf
end

function TropicalHypersurface(f::MPolyElem, val::TropicalSemiringMap, M::Union{typeof(min),typeof(max)}=min)
    tropf = tropical_polynomial(f,val)
    Tf = TropicalHypersurface(tropf)
    w = pm_object(Tf).WEIGHTS
    set_attribute!(Tf,:algebraic_polynomial,f)
    set_attribute!(Tf,:tropical_polynomial,tropf)
    set_attribute!(Tf,:weights,w)
    return Tf
end


# @doc Markdown.doc"""
#     tropical_variety(f::MPolyElem, M::Union{typeof(min),typeof(max)})

# Return the tropical variety of an algebraic polynomial in the form of a TropicalHypersurface.
# If M=min, the tropical hypersurface will obey the min-convention.
# If M=max, the tropical hypersurface will obey the max-convention.
# If coefficient ring has a valuation, the tropical hypersurface will be constructed with respect to it.
# If coefficient ring has no valuation, the tropical hypersurface will be constructed with respect to the trivial valuation.
# The function is the same as TropicalHypersurface{M}(f).

# # Examples
# julia> K = PadicField(7, 2)

# julia> Kxy, (x,y) = K["x", "y"]

# julia> f = 7*x+y+49

# julia> tropical_variety(f,min)

# julia> tropical_variety(f,max)
# """
# function tropical_variety(f::MPolyElem, M::Union{typeof(min),typeof(max)})
#     return TropicalHypersurface{M}(f)
# end

################################################################################
#
#  Basic properties for tropical hypersurfaces
#
################################################################################

# todo: add examples for varieties, curves and linear spaces
@doc Markdown.doc"""
    dual_subdivision(TH::TropicalHypersurface{M, EMB})

Return the dual subdivision of `TH` if it is embedded. Otherwise an error is thrown.

# Examples
A tropical hypersurface in $\mathbb{R}^n$ is always of dimension n-1.
```jldoctest
julia> T = TropicalSemiring(min);

julia> Txy,(x,y) = T["x","y"];

julia> f = x+y+1;

julia> tropicalLine = TropicalHypersurface(f);

julia> dual_subdivision(tropicalLine)
A subdivision of points in ambient dimension 3
```
"""
function dual_subdivision(TH::TropicalHypersurface{M,EMB}) where {M,EMB}
    if !EMB
        error("tropical hypersurface not embedded")
    end

    return SubdivisionOfPoints(pm_object(TH).DUAL_SUBDIVISION)
end
export dual_subdivision


@doc Markdown.doc"""
    polynomial(TH::TropicalHypersurface{M, EMB})

Return the tropical polynomial of `TH` if it is embedded. Otherwise an error is thrown.

# Examples
```jldoctest
julia> T = TropicalSemiring(min);

julia> Txy,(x,y) = T["x","y"];

julia> f = x+y+1;

julia> TH = TropicalHypersurface(f);

julia> polynomial(TH)
x + y + (1)
```
"""
function polynomial(TH::TropicalHypersurface{M,EMB}) where {M,EMB}
    if !EMB
        error("tropical hypersurface not embedded")
    end
    return get_attribute(TH,:tropical_polynomial)
end
export polynomial


@doc Markdown.doc"""
    minimalPolynomialFromHypersurface(T::TropicalHypersurface)

Return the minimal polynomial with smallest possible coefficients of a hypersurface.
"""
function minimalPolynomialFromHypersurface(T::TropicalHypersurface)
    return
end
