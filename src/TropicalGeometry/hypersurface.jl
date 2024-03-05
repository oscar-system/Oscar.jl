###############################################################################
#
#  Tropical hypersurfaces
#  ======================
#  concrete subtype of TropicalVarietySupertype in variety_supertype.jl
#
###############################################################################

@attributes mutable struct TropicalHypersurface{minOrMax,isEmbedded} <: TropicalVarietySupertype{minOrMax,isEmbedded}
    polyhedralComplex::PolyhedralComplex
    multiplicities::Vector{ZZRingElem}

    # tropical hypersurfaces need to be embedded
    function TropicalHypersurface{minOrMax,true}(Sigma::PolyhedralComplex, multiplicities::Vector{ZZRingElem}) where {minOrMax<:Union{typeof(min),typeof(max)}}
        @req codim(Sigma)==1 "input polyhedral complex not one-codimensional"
        return new{minOrMax,true}(Sigma,multiplicities)
    end
end



################################################################################
#
#  Printing
#
################################################################################

function Base.show(io::IO, th::TropicalHypersurface{typeof(min),true})
    print(io, "Min tropical hypersurface")
end
function Base.show(io::IO, th::TropicalHypersurface{typeof(max),true})
    print(io, "Max tropical hypersurface")
end



################################################################################
#
#  Constructors
#
################################################################################

function tropical_hypersurface(Sigma::PolyhedralComplex, mult::Vector{ZZRingElem}, minOrMax::Union{typeof(min),typeof(max)}=min)
    return TropicalHypersurface{typeof(minOrMax),true}(Sigma,mult)
end


function tropical_hypersurface(TropV::TropicalVarietySupertype{minOrMax,true}) where {minOrMax<:Union{typeof(max), typeof(min)}}
    @req codim(TropV)==1 "tropical variety codimension not one"
    @req is_pure(TropV) "tropical variety not pure"
    return tropical_hypersurface(polyhedral_complex(TropV),multiplicities(TropV),convention(TropV))
end

# Decompose and homogenize a tropical polynomial into parts that Polymake can eat.
function homogenize_and_convert_to_pm(f::Oscar.MPolyRingElem{TropicalSemiringElem{minOrMax}}) where {minOrMax<:Union{typeof(max), typeof(min)}}
   td = total_degree(f)
   exps = matrix(ZZ, collect([td-sum(alpha); alpha] for alpha in exponents(f)))
   coeffs = collect(coefficients(f))
   return coeffs, exps
end

@doc raw"""
    tropical_hypersurface(f::MPolyRingElem{<:TropicalSemiringElem}, weighted_polyhedral_complex_only::Bool=false)

Return the tropical hypersurface of the tropical polynomial `f`.  If `weighted_polyhedral_complex==true`, will not cache any extra information.

# Examples
```jldoctest
julia> T = tropical_semiring()
Min tropical semiring

julia> R,(x,y) = T["x","y"];

julia> f = x+y+1
x + y + (1)

julia> tropical_hypersurface(f)
Min tropical hypersurface

```
"""
function tropical_hypersurface(f::MPolyRingElem{<:TropicalSemiringElem}; weighted_polyhedral_complex_only::Bool=false)
    @req total_degree(f)>0 "polynomial needs to be non-constant"

    # Construct hypersurface in polymake
    minOrMax = convention(f)
    coeffs, exps = homogenize_and_convert_to_pm(f)
    pmhypproj = Polymake.tropical.Hypersurface{minOrMax}(MONOMIALS=exps, COEFFICIENTS=coeffs)
    pmhyp = Polymake.tropical.affine_chart(pmhypproj)

    # Convert to Oscar objects
    polyhedralComplex = polyhedral_complex(pmhyp)
    multiplicities = Vector{ZZRingElem}(pmhypproj.WEIGHTS)

    TropH = tropical_hypersurface(polyhedralComplex,multiplicities,minOrMax)
    if !weighted_polyhedral_complex_only
        set_attribute!(TropH,:polymake_bigobject,pmhypproj)
        set_attribute!(TropH,:tropical_polynomial,f)
    end
    return TropH
end


@doc raw"""
    tropical_hypersurface(f::MPolyRingElem, val::TropicalSemiringMap; weighted_polyhedral_complex_only::Bool=false)

Return the tropical hypersurface of the tropical polynomial that is the image of `f` under coefficient-wise `val`.  If `weighted_polyhedral_complex==true`, will not cache any extra information.

# Examples
```jldoctest
julia> R,(x,y) = QQ["x","y"];

julia> val = tropical_semiring_map(QQ,2)
Map into Min tropical semiring encoding the 2-adic valuation on Rational field

julia> f = x+y+2
x + y + 2

julia> tropical_hypersurface(f,val)
Min tropical hypersurface

```
"""
function tropical_hypersurface(f::MPolyRingElem, nu::Union{Nothing,TropicalSemiringMap}=nothing;
                               weighted_polyhedral_complex_only::Bool=false)
    # initialize nu as the trivial valuation if not specified by user
    isnothing(nu) && (nu=tropical_semiring_map(coefficient_ring(f)))

    tropf = tropical_polynomial(f,nu)
    TropH = tropical_hypersurface(tropf,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)

    if !weighted_polyhedral_complex_only
        set_attribute!(TropH,:algebraic_polynomial,f)
        set_attribute!(TropH,:tropical_semiring_map,nu)
    end
    return TropH
end


@doc raw"""
    tropical_hypersurface(Delta::SubdivisionOfPoints, minOrMax::Union{typeof(min),typeof(max)}=min; weighted_polyhedral_complex_only::Bool=false)

Construct the tropical hypersurface dual to a regular subdivision `Delta` in convention `minOrMax` using the minimal weights that give rise to it.  If `weighted_polyhedral_complex==true`, will not cache any extra information.

!!! warning
    There is a known bug when the subdivision is too easy, e.g., see example below.

# Examples
```jldoctest
julia> Delta = subdivision_of_points(simplex(2),[0,0,1])
Subdivision of points in ambient dimension 2

julia> # tropical_hypersurface(Delta) # issue 2628
```
"""
function tropical_hypersurface(Delta::SubdivisionOfPoints, minOrMax::Union{typeof(min),typeof(max)}=min;
                               weighted_polyhedral_complex_only::Bool=false)

    coeffs = min_weights(Delta)
    exps = points(Delta)
    pmhypproj = Polymake.tropical.Hypersurface{minOrMax}(MONOMIALS=exps, COEFFICIENTS=coeffs)
    pmhyp = Polymake.tropical.affine_chart(pmhypproj)

    # Convert to Oscar objects
    polyhedralComplex = polyhedral_complex(pmhyp)
    multiplicities = Vector{ZZRingElem}(pmhypproj.WEIGHTS)

    TropH = tropical_hypersurface(polyhedralComplex,multiplicities,minOrMax)
    if !weighted_polyhedral_complex_only
        set_attribute!(TropH,:polymake_object,pmhypproj)
        set_attribute!(TropH,:subdivision_of_points,Delta)
        set_attribute!(TropH,:convention,minOrMax)
    end
    return TropH
end



################################################################################
#
#  Properties
#
################################################################################
@doc raw"""
    algebraic_polynomial(TropH::TropicalHypersurface)

Return the polynomial over a valued field used to construct `TropH`.  Raises an error, if it is not cached.
"""
function algebraic_polynomial(TropH::TropicalHypersurface)
    @req has_attribute(TropH,:algebraic_polynomial) "no algebraic polynomial cached"
    return get_attribute(TropH,:algebraic_polynomial)
end


# @doc raw"""
#     convention(TropH::TropicalHypersurface)

# Return min or max depending on the convention of `TropH`.  Raises an error, if neither it nor a tropical semiring map is cached.
# """
# function convention(TropH::TropicalHypersurface)
#     if has_attribute(TropH,:convention)
#         return get_attribute(TropH,:convention)
#     elseif has_attribute(TropH,:tropical_semiring_map)
#         return convention(get_attribute(TropH,:tropical_semiring_map))
#     end
#     error("neither convention nor tropical semiring map cached")
# end


@doc raw"""
    dual_subdivision(TropH::TropicalHypersurface)

Return the dual subdivision of `TropH`.  Raises an error, if neither it nor the the internal polymake object is cached.
"""
function dual_subdivision(TropH::TropicalHypersurface{minOrMax,true}) where minOrMax
    if has_attribute(TropH,:dual_subdivision)
        return get_attribute(TropH,:dual_subdivision)
    elseif has_attribute(TropH,:polymake_bigobject)
        return get_attribute(TropH,:polymake_bigobject)
    end
    error("neither dual subdivision nor polymake object cached")
end


@doc raw"""
    tropical_polynomial(TropH::TropicalHypersurface)

Return the tropical polynomial used to construct `TropH`.  Raises an error, if it is not cached.
"""
function tropical_polynomial(TropH::TropicalHypersurface)
    @req has_attribute(TropH,:tropical_polynomial) "no tropical polynomial cached"
    return get_attribute(TropH,:tropical_polynomial)
end


# @doc raw"""
#     tropical_semiring_map(TropH::TropicalHypersurface)

# Return the tropical semiring map used to construct `TropH`.  Raises an error, if it is not cached.
# """
# function tropical_semiring_map(TropH::TropicalHypersurface)
#     @req has_attribute(TropH,:tropical_semiring_map) "no tropical semiring map cached"
#     return get_attribute(TropH,:tropical_semiring_map)
# end
