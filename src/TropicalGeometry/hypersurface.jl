###
# Tropical hypersurfaces in Oscar
###

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
        print(io, "$(repr(M)) tropical hypersurface embedded in $(ambient_dim(th))-dimensional Euclidean space")
    else
        print(io, "Abstract $(repr(M)) tropical hypersurface of dimension $(dim(th))")
    end
end

################################################################################
#
#  Constructors for tropical hypersurfaces
#
################################################################################

@doc raw"""
    TropicalHypersurface(f::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{T}})

Return the tropical hypersurface of a tropical polynomial `f`.

# Examples
```jldoctest
julia> T = TropicalSemiring(min)
Tropical semiring (min)

julia> Txy,(x,y) = T["x","y"]
(Multivariate polynomial ring in 2 variables over tropical semiring (min), AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}}[x, y])

julia> f = x+y+1
x + y + (1)

julia> Tf = TropicalHypersurface(f)
min tropical hypersurface embedded in 2-dimensional Euclidean space
```
"""
function TropicalHypersurface(f::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{T}}) where T
    if total_degree(f) <= 0
        error("Tropical hypersurfaces of constant polynomials not supported.")
    end
    M = convention(base_ring(f))
    coeffs, exps = homogenize_and_convert_to_pm(f)
    pmhypproj = Polymake.tropical.Hypersurface{M}(MONOMIALS=exps, COEFFICIENTS=coeffs)
    pmhyp = Polymake.tropical.affine_chart(pmhypproj)
    Vf = TropicalHypersurface{M, true}(polyhedral_complex(pmhyp))
    w = pmhypproj.WEIGHTS
    set_attribute!(Vf,:polymake_bigobject,pmhypproj)
    set_attribute!(Vf,:tropical_polynomial,f)
    set_attribute!(Vf,:weights,w)
    return Vf
end

# @doc raw"""
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
# A min tropical hypersurface embedded in 2-dimensional Euclidean space
# ```
# """
# function tropical_variety(f::Union{AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}},
#                                    AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(max)}}})
#     return TropicalHypersurface(f)
# end

@doc raw"""
    TropicalHypersurface(f::MPolyRingElem,M::Union{typeof(min),typeof(max)}=min)

Given a polynomial `f` over a field with an intrinsic valuation (i.e., a field
on which a function `valuation` is defined such as `PadicField(7,2)`),
return the tropical hypersurface of `f` under the convention specified by `M`.

# Examples
```jldoctest
julia> K = PadicField(7, 2);

julia> Kxy, (x,y) = K["x", "y"]
(Multivariate polynomial ring in 2 variables over QQ_7, AbstractAlgebra.Generic.MPoly{padic}[x, y])

julia> f = 7*x+y+49;

julia> TropicalHypersurface(f, min)
min tropical hypersurface embedded in 2-dimensional Euclidean space

julia> TropicalHypersurface(f, max)
max tropical hypersurface embedded in 2-dimensional Euclidean space
```
"""
function TropicalHypersurface(f::MPolyRingElem,M::Union{typeof(min),typeof(max)}=min)
    tropf = tropical_polynomial(f,M)
    Tf = TropicalHypersurface(tropf)
    w = pm_object(Tf).WEIGHTS
    set_attribute!(Tf,:algebraic_polynomial,f)
    set_attribute!(Tf,:tropical_polynomial,tropf)
    set_attribute!(Tf,:weights,w)
    return Tf
end

@doc raw"""
    TropicalHypersurface(f::MPolyRingElem,M::Union{typeof(min),typeof(max)}=min)

Construct the tropical hypersurface from a polynomial `f` and a map to the
tropical semiring `val`.

# Examples
```jldoctest
julia> Kx, (x1,x2) = polynomial_ring(QQ,2)
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x1, x2])

julia> val = TropicalSemiringMap(QQ,7)
The 7-adic valuation on Rational field

julia> f = 7*x1+x2+49;

julia> TropicalHypersurface(f, val)
min tropical hypersurface embedded in 2-dimensional Euclidean space
```
"""
function TropicalHypersurface(f::MPolyRingElem, val::TropicalSemiringMap)
    tropf = tropical_polynomial(f,val)
    Tf = TropicalHypersurface(tropf)
    w = pm_object(Tf).WEIGHTS
    set_attribute!(Tf,:algebraic_polynomial,f)
    set_attribute!(Tf,:tropical_polynomial,tropf)
    set_attribute!(Tf,:weights,w)
    return Tf
end


# @doc raw"""
#     tropical_variety(f::MPolyRingElem, M::Union{typeof(min),typeof(max)})

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
# function tropical_variety(f::MPolyRingElem, M::Union{typeof(min),typeof(max)})
#     return TropicalHypersurface{M}(f)
# end

################################################################################
#
#  Basic properties for tropical hypersurfaces
#
################################################################################

# todo: add examples for varieties, curves and linear spaces
@doc raw"""
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
Subdivision of points in ambient dimension 3
```
"""
function dual_subdivision(TH::TropicalHypersurface{M,EMB}) where {M,EMB}
    if !EMB
        error("tropical hypersurface not embedded")
    end

    return subdivision_of_points(pm_object(TH).DUAL_SUBDIVISION)
end


@doc raw"""
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


@doc raw"""
    minpoly(T::TropicalHypersurface)

Return the minimal polynomial with smallest possible coefficients of a hypersurface.
"""
function minpoly(T::TropicalHypersurface)
    error("function not implemented yet")
    return
end
