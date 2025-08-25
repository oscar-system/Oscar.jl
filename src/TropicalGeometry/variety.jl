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
    tropical_variety(I::MPolyIdeal[, nu::TropicalSemiringMap]; weighted_polyhedral_complex_only::Bool=false, skip_saturation::Bool=false, skip_decomposition::Bool=false)

Return the tropicalization of `I` with respect to `nu`.  If `nu==nothing`, will compute with respect
to the trivial valuation and min convention.  If `weighted_polyhedral_complex_only==false`, will
cache any additional information.  If `skip_saturation==false`, will saturate `I` at the product of
all variables before computing tropicalizations.  If ``skip_decomposition==false``, will return a
vector of tropical varieties, one for each primary factor of `I`.

!!! warning

    Experimental feature, only special cases supported:
    - any coefficient field and any valuation: `I` principal, binomial, or affine linear
    - QQ and trivial / p-adic valuation only: `I` primary

    Default choices for `skip_saturation` and `skip_decomposition` will change in the future to
    ensure consistency with other OSCAR functions and tropicalization functions in other software.


# Examples
```jldoctest
julia> K,t = rational_function_field(GF(101),:t);

julia> nu = tropical_semiring_map(K,t);

julia> R,(x,y,z) = K["x","y","z"];

julia> I = intersect(ideal([x+y+z+1,2*x+11*y+23*z+31]),ideal([t^3*x*y*z-1]));

julia> TropVs = tropical_variety(I,nu)
2-element Vector{TropicalVariety{typeof(min), true}}:
 Min tropical variety
 Min tropical variety

julia> K,t = rational_function_field(GF(101),:t);

julia> nu = tropical_semiring_map(K,t);

julia> R,(x,y,z) = K["x","y","z"];

julia> I = intersect(ideal([x+y+z+1,2*x+11*y+23*z+31]),ideal([t^3*x*y*z-1]));

julia> TropVs = tropical_variety(I,nu)
2-element Vector{TropicalVariety{typeof(min), true}}:
 Min tropical variety
 Min tropical variety

julia> nu_2 = tropical_semiring_map(QQ,2)
Map into Min tropical semiring encoding the 2-adic valuation on Rational field

julia> nu_3 = tropical_semiring_map(QQ,3)
Map into Min tropical semiring encoding the 3-adic valuation on Rational field

julia> f1 = 8*x^2 + x*y + x*z + x + 8*y^2 + y*z + y + 8*z^2 + z + 8;

julia> f2 = x + 2;

julia> I = ideal([f1,f2]);

julia> TropI_2 = tropical_variety(I,nu_2; skip_saturation=true, skip_decomposition=true)
Min tropical variety

julia> vertices(TropI_2)
2-element SubObjectIterator{PointVector{QQFieldElem}}:
 [-4, -4, -4]
 [-4, 0, 0]

julia> TropI_3 = tropical_variety(I,nu_3; skip_saturation=true, skip_decomposition=true)
Min tropical variety

julia> vertices(TropI_3)
1-element SubObjectIterator{PointVector{QQFieldElem}}:
 [0, -2, -2]

```
"""
function tropical_variety(I::MPolyIdeal, nu::TropicalSemiringMap=tropical_semiring_map(coefficient_ring(I)); weighted_polyhedral_complex_only::Bool=false, skip_saturation::Bool=false, skip_decomposition::Bool=false)

    if !skip_saturation
        ###
        # If saturation requested, saturate `I`
        ###
        R = base_ring(I)
        I = saturation(I,ideal([prod(gens(R))]))
    end
    # assert that `I` is not the whole ring
    @req !isone(I) "ideal contains a monomial, tropical varieties in OSCAR cannot be empty"

    if skip_decomposition
        ###
        # If decomposition not requested, return tropicalization of `I`.
        ###
        return tropical_variety_dispatch(I,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
    else
        ###
        # If decomposition requested, return tropicalization of each primary factor of `I`,
        ###
        return [tropical_variety_dispatch(P,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only) for (P,_) in primary_decomposition(I)]
    end
end


function tropical_variety(f::MPolyRingElem, nu::TropicalSemiringMap=tropical_semiring_map(coefficient_ring(f)); weighted_polyhedral_complex_only::Bool=false)
    return tropical_variety(ideal(parent(f),[f]),nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
end

###
# Function for dispatching which tropicalization routine to call
###
function tropical_variety_dispatch(I::MPolyIdeal, nu::TropicalSemiringMap=tropical_semiring_map(coefficient_ring(I)); weighted_polyhedral_complex_only::Bool=false)
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

    # I general
    return tropical_variety_prime(I,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
end
