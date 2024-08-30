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

function Base.show(io::IO, tv::TropicalVariety{typeof(min), true})
    print(io, "Min tropical variety")
end
function Base.show(io::IO, tv::TropicalVariety{typeof(max), true})
    print(io, "Max tropical variety")
end



################################################################################
#
#  Constructors
#
################################################################################

@doc raw"""
    tropical_variety(Sigma::PolyhedralComplex, mult, minOrMax::Union{typeof(min),typeof(max)}=min)

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

    # copy all attributes
    TropV.__attrs = deepcopy(TropH.__attrs)

    return TropV
end


function tropical_variety(TropL::TropicalLinearSpace{minOrMax,true}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    TropV = tropical_variety(polyhedral_complex(TropL),multiplicities(TropL), convention(TropL))

    # copy all attributes
    TropV.__attrs = deepcopy(TropL.__attrs)

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
    tropical_variety(I::MPolyIdeal, nu::Union{TropicalSemiringMap,Nothing}=nothing; weighted_polyhedral_complex_only::Bool=false)

Return the tropicalization of `I` with respect to `nu`.
If `nu==nothing`, will compute with respect to the trivial valuation and min convention.
If `weighted_polyhedral_complex_only==true`, will not cache any additional information.

!!! warning
    Assumes that `I` is equi-dimensional.  Only special cases supported:
    - any valuation: `I` principal, binomial, affine linear
    - trivial and p-adic valuation only: `I` general

# Examples
```jldoctest
julia> R,(x,y) = QQ[:x, :y];

julia> I = ideal([(x^2+y)*(x+y^2)*(x+y)]);

julia> tropical_variety(I)
3-element Vector{TropicalVariety}:
 Min tropical variety
 Min tropical variety
 Min tropical variety

```
"""
function tropical_variety(I::MPolyIdeal, nu::Union{TropicalSemiringMap,Nothing}=nothing; weighted_polyhedral_complex_only::Bool=false)
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

    # I general
    return tropical_variety_equidimensional(I,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
end


function tropical_variety(f::MPolyRingElem, nu::Union{TropicalSemiringMap,Nothing}=nothing; weighted_polyhedral_complex_only::Bool=false)
    return tropical_variety(ideal(parent(f),[f]),nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
end


function homogenize_pre_tropicalization(I::MPolyIdeal)
    ###
    # Compute reduced Groebner basis (usually already cached), and construct homogenization
    ###
    G = groebner_basis(I,complete_reduction=true)

    Kx = base_ring(I)
    K = coefficient_ring(Kx)
    x = symbols(Kx)
    Kxhx,_ = polynomial_ring(K,vcat([:xh],x))

    Gh = Vector{elem_type(Kx)}(undef,length(G))
    for (i,g) in enumerate(G)
        gh = MPolyBuildCtx(Kxhx)
        d = total_degree(g)
        for (c,alpha) in coefficients_and_exponents(g)
            pushfirst!(alpha,d-sum(alpha)) # homogenize exponent vector
            push_term!(gh,c,alpha)
        end
        Gh[i] = finish(gh)
    end

    return ideal(Gh)
end


function dehomogenize_post_tropicalization(Sigma::PolyhedralComplex)
    @req lineality_dim(Sigma)>0 "dehomogenizing polyhedral complex without lineality"

    ###
    # Construct hyperplane {first coord = 0}
    ###
    n = ambient_dim(Sigma)
    zerothUnitRowVector = zeros(Int,1,n)
    zerothUnitRowVector[1,1] = 1
    dehomogenisingHyperplane = polyhedron((zeros(Int,0,n),zeros(Int,0)), (zerothUnitRowVector,[0]))

    ###
    # Construct matrix and incidence matrix of vertices and rays
    ###
    incidenceMatrixVertices = Vector{Int}[]
    dehomogenizedVertices = Vector{QQFieldElem}[]
    incidenceMatrixRays = Vector{Int}[]
    dehomogenizedRays = Vector{QQFieldElem}[]
    for sigma in maximal_polyhedra(Sigma)
        sigmaDehomogenized = intersect(sigma,dehomogenisingHyperplane)
        incidenceVectorVertices = Int[]
        V,_ = minimal_faces(sigmaDehomogenized)
        for vertex in V
            vertex = vertex[2:end]
            i = findfirst(isequal(vertex),dehomogenizedVertices)
            if i === nothing
                push!(dehomogenizedVertices,vertex)
                push!(incidenceVectorVertices,length(dehomogenizedVertices))
            else
                push!(incidenceVectorVertices,i)
            end
        end
        push!(incidenceMatrixVertices,incidenceVectorVertices)

        incidenceVectorRays = Int[]
        R,_ = rays_modulo_lineality(sigmaDehomogenized)
        for ray in R
            ray = ray[2:end]
            i = findfirst(isequal(ray),dehomogenizedRays)
            if i === nothing
                push!(dehomogenizedRays,ray)
                push!(incidenceVectorRays,length(dehomogenizedRays))
            else
                push!(incidenceVectorRays,i)
            end
        end
        push!(incidenceMatrixRays,incidenceVectorRays)
    end

    ###
    # Concatenate vertically matrixes of vertices and rays,
    # shift incidence matrix of rays and concatenate it horizontally to incicende matrix of vertices,
    # dehomogenize generators of lineality space
    ###
    dehomogenizedVerticesAndRays = matrix(QQ,vcat(dehomogenizedVertices,dehomogenizedRays))
    incidenceMatrixRaysShifted = (x -> x .+length(dehomogenizedVertices)).(incidenceMatrixRays)
    incidenceMatrixVerticesAndRays = IncidenceMatrix([vcat(iv,ir) for (iv,ir) in zip(incidenceMatrixVertices,incidenceMatrixRaysShifted)])

    ###
    # Dehomogenize lineality space
    ###
    sigma = first(maximal_polyhedra(Sigma))
    sigmaDehomogenized = intersect(sigma,dehomogenisingHyperplane)
    dehomogenizedLineality = [linealityVector[2:end] for linealityVector in lineality_space(sigmaDehomogenized)]

    return polyhedral_complex(incidenceMatrixVerticesAndRays,
                              dehomogenizedVerticesAndRays,
                              collect(length(dehomogenizedVertices)+1:length(dehomogenizedVertices)+length(dehomogenizedRays)),
                              dehomogenizedLineality)
end
