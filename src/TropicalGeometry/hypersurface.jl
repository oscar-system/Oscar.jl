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
    @req codim(TropV)==1 "tropical variety not one-codimensional"
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

julia> R,(x,y) = T[:x, :y];

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
julia> R,(x,y) = QQ[:x, :y];

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

Construct the tropical hypersurface dual to a regular subdivision `Delta` in convention `minOrMax`. To be precise, the tropical hypersurface of the tropical polynomial with exponent vectors `points(Delta)` and coefficients `min_weight(Delta)` (min-convention) or `-min_weight(Delta)` (max-convention).  If `weighted_polyhedral_complex==true`, will not cache any extra information.

# Examples
```jldoctest
julia> Delta = subdivision_of_points([0 0; 1 0; 0 1; 2 0],[0,0,0,1])
Subdivision of points in ambient dimension 2

julia> tropical_hypersurface(Delta)
Min tropical hypersurface
```
"""
function tropical_hypersurface(Delta::SubdivisionOfPoints, minOrMax::Union{typeof(min),typeof(max)}=min;
                               weighted_polyhedral_complex_only::Bool=false)

    coeffs = min_weights(Delta)
    expvs = [ZZ.(alpha) for alpha in points(Delta)]

    TT = tropical_semiring(minOrMax)
    # preserve_ordering=true, since min_weights of regular subdivisions are always in min-convention,
    # e.g., [0 0; 1 0; 0 1; 2 0] decomposed into [0 0; 1 0; 0 1] and [1 0; 0 1; 2 0] has min_weight [+1,0,0,0]
    # which is dual to the tropical hypersurface of min(+1, x, y, 2*x) or max(-1, x, y, 2*x)
    coeffs = TT.(coeffs; preserve_ordering=true)
    _,x = polynomial_ring(TT,length(first(expvs)); cached=false)
    tropf = sum([c*prod(x.^alpha) for (c,alpha) in zip(coeffs,expvs)])
    TropH = tropical_hypersurface(tropf)

    if !weighted_polyhedral_complex_only
        set_attribute!(TropH,:dual_subdivision,Delta)
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

Return the polynomial over a valued field used to construct `TropH`.  Raises an error if it is not cached.
"""
function algebraic_polynomial(TropH::TropicalHypersurface)
    @req has_attribute(TropH,:algebraic_polynomial) "no algebraic polynomial cached"
    return get_attribute(TropH, :algebraic_polynomial)::MPolyRingElem
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

Return the dual subdivision used to construct `TropH`.  Raises an error if it is not cached.

# Examples
```jldoctest
julia> Delta = subdivision_of_points([0 0; 1 0; 0 1; 2 0],[0,0,0,1])
Subdivision of points in ambient dimension 2

julia> th = tropical_hypersurface(Delta)
Min tropical hypersurface

julia> sop = dual_subdivision(th)
Subdivision of points in ambient dimension 2

julia> points(sop)
4-element SubObjectIterator{PointVector{QQFieldElem}}:
 [0, 0]
 [1, 0]
 [0, 1]
 [2, 0]

julia> maximal_cells(sop)
2-element SubObjectIterator{Vector{Int64}}:
 [1, 2, 3]
 [2, 3, 4]
```
"""
function dual_subdivision(TropH::TropicalHypersurface{minOrMax,true}) where minOrMax
    @req has_attribute(TropH,:dual_subdivision) "no dual subdivision cached"
    return get_attribute(TropH, :dual_subdivision)::SubdivisionOfPoints
end


@doc raw"""
    tropical_polynomial(TropH::TropicalHypersurface)

Return the tropical polynomial used to construct `TropH`.  Raises an error if it is not cached.
"""
function tropical_polynomial(TropH::TropicalHypersurface{minOrMax}) where minOrMax
    @req has_attribute(TropH,:tropical_polynomial) "no tropical polynomial cached"
    return get_attribute(TropH,:tropical_polynomial)::mpoly_type(TropicalSemiringElem{minOrMax})
end


# @doc raw"""
#     tropical_semiring_map(TropH::TropicalHypersurface)

# Return the tropical semiring map used to construct `TropH`.  Raises an error, if it is not cached.
# """
# function tropical_semiring_map(TropH::TropicalHypersurface)
#     @req has_attribute(TropH,:tropical_semiring_map) "no tropical semiring map cached"
#     return get_attribute(TropH,:tropical_semiring_map)
# end
#

@doc raw"""
   tropical_prevariety(F::Vector{MPolyRingElem},nu::TropicalSemiringMap)

Return the tropical prevariety generated by intersecting tropical hypersurfaces corresponding to elements of `F`.

If `F` is a finite collection of polynomials with coefficients from a given field, return the tropical prevariety obtained by tropicalizing elements of `F` with respect to a given tropicalization map `nu`.

If no `nu` is given, default to trivial valuation with `min` convention.

If `F` is a collection of tropical polynomials, the function computes and intersects the associated hypersurfaces.
# Example
We compute the Dressian $\text{Dr}(2,5)$ below.
```jldoctest
julia> Gr25 = grassmann_pluecker_ideal(2,5)
Ideal generated by
  x[[1, 2]]*x[[3, 4]] - x[[1, 3]]*x[[2, 4]] + x[[1, 4]]*x[[2, 3]]
  x[[1, 2]]*x[[3, 5]] - x[[1, 3]]*x[[2, 5]] + x[[1, 5]]*x[[2, 3]]
  x[[1, 2]]*x[[4, 5]] - x[[1, 4]]*x[[2, 5]] + x[[1, 5]]*x[[2, 4]]
  x[[1, 3]]*x[[4, 5]] - x[[1, 4]]*x[[3, 5]] + x[[1, 5]]*x[[3, 4]]
  x[[2, 3]]*x[[4, 5]] - x[[2, 4]]*x[[3, 5]] + x[[2, 5]]*x[[3, 4]]

#Compute Dressian without specified tropicalization map.
julia> Dr25 = tropical_prevariety(gens(Gr25))
Polyhedral complex in ambient dimension 10

julia> rays_modulo_lineality(Dr25)[1]
10-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, -1//3, -1//3, -1//3, -1//3, -1//3, -1//3, 1//3, 1//3, 1//3]
 [1, 1, -1, -1, 1, -1, -1, -1, -1, 3]
 [1, -1, 1, -1, -1, 1, -1, -1, 3, -1]
 [1, -1, -1, 1, -1, -1, 1, 3, -1, -1]
 [-1, 1, 1, -1, -1, -1, 3, 1, -1, -1]
 [-1, -1, -1, 3, 1, 1, -1, 1, -1, -1]
 [-1, 1, -1, 1, -1, 3, -1, -1, 1, -1]
 [-1, 3, -1, -1, -1, 1, 1, -1, -1, 1]
 [-1, -1, 3, -1, 1, -1, 1, -1, 1, -1]
 [-1, -1, 1, 1, 3, -1, -1, -1, -1, 1]

#Compute with respect to max convention
julia> nu = tropical_semiring_map(QQ,max)
Map into Max tropical semiring encoding the trivial valuation on Rational field

julia> Dr25max = tropical_prevariety(gens(Gr25),nu)
Polyhedral complex in ambient dimension 10

#Give input as tropical polynomials
Dr25 = tropical_prevariety(tropical_polynomial.([f for f in gens(Gr25)]))
Polyhedral complex in ambient dimension 10
```
"""
#Input array of tropical polynomials
function tropical_prevariety(F::Vector{<:AbstractAlgebra.Generic.MPoly{TropicalSemiringElem{minmax}}}) where {minmax<:Union{typeof(min),typeof(max)}}
    length(F) > 0 || error("Empty array.")
    HTF = tropical_hypersurface.(F)
    return tropical_prevariety(HTF)
end

function tropical_prevariety(F::Vector{<:MPolyRingElem}, nu::TropicalSemiringMap)
   return tropical_prevariety(tropical_polynomial.(F,Ref(nu)))
end

#Default to min with trivial valuation
tropical_prevariety(F::Vector{<:MPolyRingElem}) = tropical_prevariety(F, tropical_semiring_map(coefficient_ring(first(F))))

tropical_prevariety(F::Vector{<:MPolyRingElem}, convention::minmax) where {minmax<:Union{typeof(min),typeof(max)}} = tropical_prevariety(F, tropical_semiring_map(coefficient_ring(first(F)), convention))

tropical_prevariety(F::Vector{TropicalHypersurface{minmax, true}}) where {minmax<:Union{typeof(min),typeof(max)}} = reduce(common_refinement, [h.polyhedralComplex for h in F])
