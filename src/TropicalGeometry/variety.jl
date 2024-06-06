################################################################################
#
#  Tropical varieties
#  ==================
#  concrete subtype of TropicalVarietySupertype in variety_supertype.jl
#
################################################################################

@attributes mutable struct TropicalVariety{minOrMax,isEmbedded} <: TropicalVarietySupertype{minOrMax,isEmbedded}
    polyhedralComplex::PolyhedralComplex
    multiplicities::Vector{ZZRingElem}

    # tropical varieties need to be embedded
    function TropicalVariety{minOrMax,true}(Sigma::PolyhedralComplex,multiplicities::Vector{ZZRingElem}) where {minOrMax<:Union{typeof(min),typeof(max)}}
        return new{minOrMax,true}(Sigma,multiplicities)
    end
end



################################################################################
#
#  Printing
#
################################################################################

function Base.show(io::IO, ::TropicalVariety{typeof(min), true})
    print(io, "Min tropical variety")
end
function Base.show(io::IO, ::TropicalVariety{typeof(max), true})
    print(io, "Max tropical variety")
end



################################################################################
#
#  Constructors
#
################################################################################

function copy_tropical_attributes!(TropVcopyTo::TropicalVarietySupertype, TropVcopyFrom::TropicalVarietySupertype)
    inheritableAttributes =
        [
            # from all
            :tropical_semiring_map,
            # from tropical hypersurfaces
            :tropical_polynomial,
            :algebraic_polynomial,
            :dual_subdivision,
            # from tropical linear spaces
            :polymake_object,
            :pluecker_indices,
            :tropical_pluecker_vector,
            :tropical_matrix,
            :algebraic_pluecker_vector,
            :algebraic_matrix,
            :algebraic_ideal
            # from tropical varieties
        ]
    for attr in inheritableAttributes
        if has_attribute(TropVcopyFrom,attr)
            set_attribute!(TropVcopyTo,attr,get_attribute(TropVcopyFrom,attr))
        end
    end
end

@doc raw"""
    tropical_variety(Sigma::PolyhedralComplex, mult::Vector{ZZRingElem}, minOrMax::Union{typeof(min),typeof(max)}=min)

Return the `TropicalVariety` whose polyhedral complex is `Sigma` with multiplicities `mult` and convention `minOrMax`. Here, `mult` is optional can be specified as a `Vector{ZZRingElem}` which represents a list of multiplicities on the maximal polyhedra in the order of `maximal_polyhedra(Sigma)`.  If `mult` is unspecified, then all multiplicities are set to one.

# Examples
```jldoctest
julia> Sigma = polyhedral_complex(incidence_matrix([[1],[2]]), [[0],[1]])
Polyhedral complex in ambient dimension 1

julia> tropical_variety(Sigma)
Min tropical variety

julia> mult = ones(ZZRingElem, n_maximal_polyhedra(Sigma))
2-element Vector{ZZRingElem}:
 1
 1

julia> tropical_variety(Sigma,mult,min)
Min tropical variety

julia> mult = ZZ.([1,2])
2-element Vector{ZZRingElem}:
 1
 2

julia> tropical_variety(Sigma,mult,max)
Max tropical variety

```
"""
function tropical_variety(Sigma::PolyhedralComplex, mult::Vector{ZZRingElem}, minOrMax::Union{typeof(min),typeof(max)}=min)
    return TropicalVariety{typeof(minOrMax),true}(Sigma,mult)
end
function tropical_variety(Sigma::PolyhedralComplex, minOrMax::Union{typeof(min),typeof(max)}=min)
    mult = ones(ZZRingElem, n_maximal_polyhedra(Sigma))
    return tropical_variety(Sigma,mult,minOrMax)
end

function tropical_variety(TropH::TropicalHypersurface{minOrMax,true}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    TropV = tropical_variety(polyhedral_complex(TropH),multiplicities(TropH), convention(TropH))
    copy_tropical_attributes!(TropV,TropH)
    return TropV
end

function tropical_variety(TropL::TropicalLinearSpace{minOrMax,true}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    TropV = tropical_variety(polyhedral_complex(TropL),multiplicities(TropL), convention(TropL))
    copy_tropical_attributes!(TropV,TropL)
    return TropV
end


function tropical_variety(Sigma::Vector{<:Polyhedron}, multiplicities::Vector{ZZRingElem}, minOrMax::Union{typeof(min),typeof(max)}=min)
    return tropical_variety(polyhedral_complex(Sigma; non_redundant=true),multiplicities,minOrMax)
end




################################################################################
#
#  Properties
#  ----------
#  none, see variety_supertype.jl for common properties of all tropical variety types
#
################################################################################



################################################################################
#
#  Tropical varieties of polynomial ideals
#  ----------
#  References for computing tropical varieties via Groebner complex traversal:
#    T. Bogart, A. Jensen, D. Speyer, B. Sturmfels, R. Thomas: Computing tropical varieties
#    T. Markwig, Y. Ren: Computing tropical points over fields with valuation
#
################################################################################

@doc raw"""
    tropical_variety(I::MPolyIdeal, nu::Union{TropicalSemiringMap,Nothing}=nothing; weighted_polyhedral_complex_only::Bool=false, check::Bool=true)

Return the tropicalization of `I` with respect to `nu`.
If `nu==nothing`, will compute with respect to the trivial valuation and min convention.
If `weighted_polyhedral_complex_only==true`, will not cache any additional information.

If `check==true` and `I` is neither principal, binomial, nor affine linear, will check whether `I` is primary.

!!! warning
    Experimental feature, only special cases supported:
    - any coefficient field and any valuation: `I` principal, binomial, or affine linear
    - QQ and trivial / p-adic valuation only: `I` primary

# Examples
```jldoctest
julia> R,(x,y,z) = QQ["x","y","z"];

julia> nu_2 = tropical_semiring_map(QQ,2)
Map into Min tropical semiring encoding the 2-adic valuation on Rational field

julia> f1 = 8*x^2 + x*y + x*z + x + 8*y^2 + y*z + y + 8*z^2 + z + 8;

julia> f2 = x + 1;

julia> I = ideal([f1,f2]);

julia> TropI_0 = tropical_variety(I)
Min tropical variety

julia> vertices(TropI_0)
1-element SubObjectIterator{PointVector{QQFieldElem}}:
 [0, 0, 0]

julia> TropI_2 = tropical_variety(I,nu_2)
Min tropical variety

julia> vertices(TropI_2)
2-element SubObjectIterator{PointVector{QQFieldElem}}:
 [0, -3, 3]
 [0, 3, -3]

```
"""
function tropical_variety(I::MPolyIdeal, nu::Union{TropicalSemiringMap,Nothing}=nothing; weighted_polyhedral_complex_only::Bool=false, check::Bool=true)

    # initialize nu as trivial valuation if not given
    if isnothing(nu)
        nu = tropical_semiring_map(coefficient_ring(I))
    end

    # compute a reduced GB to test whether `I` is principal, binomial, or affine linear
    GB = groebner_basis(I,complete_reduction=true)
    if length(GB)==1
        # I is principal
        return tropical_variety_principal(I,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
    elseif all(g->(length(g)==2), GB)
        # I is binomial
        return tropical_variety_binomial(I,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
    elseif all(g->(total_degree(g)<2), GB)
        # I affine linear
        return tropical_variety_affine_linear(I,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
    end

    @req !check || is_primary(I) "Input ideal not primary. See `tropical_varieties` for non-primary ideals or disable check."

    # I general
    return tropical_variety_prime(I,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
end


function tropical_variety(f::MPolyRingElem, nu::Union{TropicalSemiringMap,Nothing}=nothing; weighted_polyhedral_complex_only::Bool=false)
    return tropical_variety(ideal(parent(f),[f]),nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
end


@doc raw"""
    tropical_varieties(I::MPolyIdeal, nu::Union{TropicalSemiringMap,Nothing}=nothing; weighted_polyhedral_complex_only::Bool=false)

Compute the primary decomposition of `I` and the tropicalizations of each primary factor.

# Examples
```jldoctest
julia> K,t = rational_function_field(GF(101),:t);

julia> nu = tropical_semiring_map(K,t);

julia> R,(x,y,z) = K["x","y","z"];

julia> I = intersect(ideal([x+y+z+1,2*x+11*y+23*z+31]),ideal([t^3*x*y*z-1]));

julia> TropVs = tropical_varieties(I,nu)
2-element Vector{TropicalVariety{typeof(min), true}}:
 Min tropical variety
 Min tropical variety

```
"""
function tropical_varieties(I::MPolyIdeal, nu::Union{TropicalSemiringMap,Nothing}=nothing; weighted_polyhedral_complex_only::Bool=false)
    PQs = primary_decomposition(I)
    Ps = first.(PQs)
    return [tropical_variety(P,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only) for P in Ps]
end
