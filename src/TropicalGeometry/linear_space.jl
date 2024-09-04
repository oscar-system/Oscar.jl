###############################################################################
#
#  Tropical linear spaces
#  ======================
#  concrete subtype of TropicalVarietySupertype in variety_supertype.jl
#
###############################################################################

@attributes mutable struct TropicalLinearSpace{minOrMax,isEmbedded} <: TropicalVarietySupertype{minOrMax,isEmbedded}
    polyhedralComplex::PolyhedralComplex
    multiplicities::Vector{ZZRingElem}

    # tropical linear spaces need to be embedded
    # Note: multiplicities needs to be passed explicitly to avoid complications with a third optional parameter (tropical semiring map)
    function TropicalLinearSpace{minOrMax,true}(Sigma::PolyhedralComplex, multiplicities::Vector{ZZRingElem}) where {minOrMax<:Union{typeof(min),typeof(max)}}
        return new{minOrMax,true}(Sigma,multiplicities)
    end
end



###############################################################################
#
#  Printing
#
###############################################################################

function Base.show(io::IO, th::TropicalLinearSpace{typeof(min), true})
    print(io, "Min tropical linear space")
end
function Base.show(io::IO, th::TropicalLinearSpace{typeof(max), true})
    print(io, "Max tropical linear space")
end



###############################################################################
#
#  Constructors
#
###############################################################################

function tropical_linear_space(Sigma::PolyhedralComplex, mult::Vector{ZZRingElem}, minOrMax::Union{typeof(min),typeof(max)}=min)
    return TropicalLinearSpace{typeof(minOrMax),true}(Sigma,mult)
end


function tropical_linear_space(TropV::TropicalVarietySupertype{minOrMax,true}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    mult = multiplicities(TropV)
    @req isnothing(findfirst(!isequal(one(ZZ)),collect(values(mult)))) "tropical variety not all multiplicities one"
    return tropical_linear_space(polyhedral_complex(TropV),mult,convention(TropV))
end


function shift_pluecker_indices_for_polymake(plueckerIndices::Vector{Vector{Int}})
    # substract 1 fromt all entries of plueckerIndices
    return (x -> x .-1).(plueckerIndices)
end


function add_missing_lineality_from_polymake(Sigma::PolyhedralComplex)
    # Adding missing ones vector to lineality space
    IM = maximal_polyhedra(IncidenceMatrix,Sigma)
    VR = vertices_and_rays(Sigma)
    rayIndices = findall(vr->(vr isa RayVector),VR)
    VRmat = matrix(QQ,Vector{QQFieldElem}.(collect(VR)))
    L = vcat(matrix(QQ,lineality_space(Sigma)),
             matrix(QQ,ones(Int,1,ambient_dim(Sigma))))
    return polyhedral_complex(IM,VRmat,rayIndices,L; non_redundant=true)
end


@doc raw"""
    tropical_linear_space(Lambda::Vector{Vector{Int}}, p::Vector{<:TropicalSemiringElem}; weighted_polyhedral_complex_only::Bool=false)

Return a tropical linear space from a tropical Pluecker vector with indices `Lambda` and values `p`.  If `weighted_polyhedral_complex==true`, will not cache any extra information.

# Examples
```jldoctest
julia> T = tropical_semiring();

julia> plueckerIndices = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]];

julia> plueckerVector = T.([0,0,0,0,0,0]);

julia> tropical_linear_space(plueckerIndices,plueckerVector)
Min tropical linear space

```
"""
function tropical_linear_space(plueckerIndices::Vector{Vector{Int}}, plueckerVector::Vector{<:TropicalSemiringElem}; weighted_polyhedral_complex_only::Bool=false)

    # read off convention and indices of finite entries
    minOrMax = convention(plueckerVector)
    plueckerSupport = findall(!iszero,plueckerVector)

    # restrict plueckerVector and pluckerIndices,
    # convert plueckerVector to QQ (for Polymake)
    plueckerVectorForPolymake = QQ.(plueckerVector[plueckerSupport])
    plueckerIndicesForPolymake = shift_pluecker_indices_for_polymake(plueckerIndices[plueckerSupport])

    # find out ambient dimension from remaining indices
    n = length(unique(Iterators.flatten(plueckerIndices)))

    # construct tropical linear space
    valuatedMatroid = Polymake.matroid.ValuatedMatroid{minOrMax}(
        BASES = plueckerIndicesForPolymake,
        N_ELEMENTS = n,
        VALUATION_ON_BASES = plueckerVectorForPolymake)

    Sigma = PolyhedralComplex{QQFieldElem}(Polymake.tropical.linear_space{minOrMax}(valuatedMatroid))
    Sigma = add_missing_lineality_from_polymake(Sigma)
    multiplicities = ones(ZZRingElem, n_maximal_polyhedra(Sigma))

    TropL = tropical_linear_space(Sigma,multiplicities,minOrMax)
    if !weighted_polyhedral_complex_only
        set_attribute!(TropL,:polymake_object,valuatedMatroid)
        set_attribute!(TropL,:pluecker_indices,plueckerIndices)
        set_attribute!(TropL,:tropical_pluecker_vector,plueckerVector)
        set_attribute!(TropL,:valuated_matroid,valuatedMatroid)
    end
    return TropL
end


@doc raw"""
    tropical_linear_space(k::Int, n::Int, p::Vector{<:TropicalSemiringElem}; weighted_polyhedral_complex_only::Bool=false)

Return a tropical linear space from a tropical Pluecker vector with indices `AbstractAlgebra.combinations(1:n,k)` and values `p`.  If `weighted_polyhedral_complex==true`, will not cache any extra information.

# Examples
```jldoctest
julia> T = tropical_semiring();

julia> plueckerVector = T.([0,0,0,0,0,0]);

julia> tropical_linear_space(2,4,plueckerVector)
Min tropical linear space

```
"""
function tropical_linear_space(k::Int, n::Int, plueckerVector::Vector{<:TropicalSemiringElem}; weighted_polyhedral_complex_only::Bool=false)
    return tropical_linear_space(AbstractAlgebra.combinations(1:n,k), plueckerVector, weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
end


@doc raw"""
    tropical_linear_space(Lambda::Vector{Vector{Int}}, p::Vector, nu::TropicalSemiringMap; weighted_polyhedral_complex_only::Bool=false)

Return a tropical linear space from a tropical Pluecker vector with indices `Lambda` and values `nu(p)`.  If `weighted_polyhedral_complex==true`, will not cache any extra information.

# Examples
```jldoctest
julia> nu = tropical_semiring_map(QQ,2)
Map into Min tropical semiring encoding the 2-adic valuation on Rational field

julia> plueckerIndices = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]];

julia> plueckerVector = QQ.([1,3,5,5,3,1]);

julia> tropical_linear_space(plueckerIndices,plueckerVector,nu)
Min tropical linear space

```
"""
function tropical_linear_space(plueckerIndices::Vector{Vector{Int}}, plueckerVector::Vector, nu::Union{Nothing,TropicalSemiringMap}=nothing; weighted_polyhedral_complex_only::Bool=false)
    # if nu unspecified, initialize as the trivial valuation + min convention
    isnothing(nu) && (nu=tropical_semiring_map(parent(first(plueckerVector))))

    TropL = tropical_linear_space(plueckerIndices,
                                  nu.(plueckerVector),
                                  weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)

    if !weighted_polyhedral_complex_only
        set_attribute!(TropL,:algebraic_pluecker_vector,plueckerVector)
        set_attribute!(TropL,:tropical_semiring_map,nu)
    end
    return TropL
end


@doc raw"""
    tropical_linear_space(k::Int, n::Int, p::Vector, nu::TropicalSemiringMap; weighted_polyhedral_complex_only::Bool=false)

Return a tropical linear space from a tropical Pluecker vector with indices `AbstractAlgebra.combinations(1:n,k)` and values `nu(p)`.  If `weighted_polyhedral_complex==true`, will not cache any extra information.

# Examples
```jldoctest
julia> nu = tropical_semiring_map(QQ);

julia> plueckerVector = QQ.([1,3,5,5,3,1]);

julia> tropical_linear_space(2,4,plueckerVector,nu)
Min tropical linear space

```
"""
function tropical_linear_space(k::Int, n::Int, plueckerVector::Vector, nu::Union{Nothing,TropicalSemiringMap}=nothing; weighted_polyhedral_complex_only::Bool=false)
    # if nu unspecified, initialize as the trivial valuation + min convention
    isnothing(nu) && (nu=tropical_semiring_map(parent(first(plueckerVector))))

    TropL = tropical_linear_space(AbstractAlgebra.combinations(1:n,k), nu.(plueckerVector), weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)

    if !weighted_polyhedral_complex_only
        set_attribute!(TropL,:algebraic_pluecker_vector,plueckerVector)
        set_attribute!(TropL,:tropical_semiring_map,nu)
    end
    return TropL
end


# helper function that takes a matrix (tropical or algebraic)
# and returns the non-zero minors and their indices
function compute_pluecker_indices_and_vector(A::MatElem)
    n = ncols(A)
    k = nrows(A)
    @req n>=k "matrix for Pluecker vector cannot have more rows than columns"
    plueckerIndices = AbstractAlgebra.combinations(1:n,k)
    plueckerVector = AbstractAlgebra.minors(A,k)
    nonZeroIndices = findall(!iszero,plueckerVector)
    plueckerIndices = plueckerIndices[nonZeroIndices]
    plueckerVector = plueckerVector[nonZeroIndices]
    return plueckerIndices, plueckerVector
end


@doc raw"""
    tropical_linear_space(A::MatElem{<:TropicalSemiringElem}; weighted_polyhedral_complex_only::Bool=false)

Return a tropical linear space whose Pluecker vector are the tropical minors of `A`.  Assumes that `ncols(A)>=nrows(A)`.  If `weighted_polyhedral_complex==true`, will not cache any extra information.

# Examples
```jldoctest
julia> A = matrix(tropical_semiring(),[[1,2,4,8],[8,4,2,1]])
[(1)   (2)   (4)   (8)]
[(8)   (4)   (2)   (1)]

julia> tropical_linear_space(A)
Min tropical linear space

```
"""
function tropical_linear_space(A::MatElem{T}; weighted_polyhedral_complex_only::Bool=false) where {T<:TropicalSemiringElem}
    @req ncols(A)>=nrows(A) "matrix for Pluecker vector cannot have more rows than columns"

    plueckerIndices,plueckerVector = compute_pluecker_indices_and_vector(A)
    TropL = tropical_linear_space(plueckerIndices, plueckerVector,
                                  weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
    if !weighted_polyhedral_complex_only
        set_attribute!(TropL,:tropical_matrix,A)
    end
    return TropL
end
function tropical_linear_space(A::Matrix{T}; weighted_polyhedral_complex_only::Bool=false) where {T<:TropicalSemiringElem}
    # convert to Oscar matrix (necessary for AbstractAlgebra.minors)
    return tropical_linear_space(matrix(parent(first(A)),A), weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
end

@doc raw"""
    tropical_linear_space(A::MatElem,nu::TropicalSemiringMap; weighted_polyhedral_complex_only::Bool=false)

Return a tropical linear space whose Pluecker vector is `nu` applied to the minors of `A`.  If `weighted_polyhedral_complex==true`, will not cache any extra information.

# Examples
```jldoctest
julia> nu = tropical_semiring_map(QQ,2)
Map into Min tropical semiring encoding the 2-adic valuation on Rational field

julia> A = matrix(QQ,[[1,2,4,8],[8,4,2,1]])
[1   2   4   8]
[8   4   2   1]

julia> tropical_linear_space(A, nu)
Min tropical linear space

```
"""
function tropical_linear_space(A::MatElem, nu::Union{Nothing,TropicalSemiringMap}=nothing; weighted_polyhedral_complex_only::Bool=false)
    # compute reduced row echelon form of A
    # and remove all zero rows so that matrix is of full rank
    _,A = rref(A)
    nonzeroRowIndices = findall(!iszero,[A[i,:] for i in 1:nrows(A)])
    A = A[nonzeroRowIndices,:]

    plueckerIndices,plueckerVector = compute_pluecker_indices_and_vector(A)
    TropL = tropical_linear_space(plueckerIndices, plueckerVector, nu,
                                  weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
    if !weighted_polyhedral_complex_only
        set_attribute!(TropL,:algebraic_matrix,A)
        set_attribute!(TropL,:tropical_semiring_map,nu)
    end
    return TropL
end
function tropical_linear_space(A::Matrix, nu::Union{Nothing,TropicalSemiringMap}=nothing; weighted_polyhedral_complex_only::Bool=false)
    return tropical_linear_space(matrix(parent(first(A)),A), nu, weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
end

# Special implementation for trivial valuation over the rationals using Polymake.tropical.matroid_fan for better performance
function tropical_linear_space(A::QQMatrix, nu::TropicalSemiringMap{QQField, Nothing, minOrMax}=tropical_semiring_map(QQ); weighted_polyhedral_complex_only::Bool=false) where {minOrMax<:Union{typeof(min),typeof(max)}}
    # compute reduced row echelon form of A
    # and remove all zero rows so that matrix is of full rank
    _,A = rref(A)
    nonzeroRowIndices = findall(!iszero,[A[i,:] for i in 1:nrows(A)])
    A = A[nonzeroRowIndices,:]

    # Compute Bergman fan in polymake and convert to polyhedral complex
    Sigma = PolyhedralComplex{QQFieldElem}(Polymake.tropical.matroid_fan{convention(nu)}(A))

    # Add the missing all ones vector to the lineality
    Sigma = add_missing_lineality_from_polymake(Sigma)
    multiplicities = ones(ZZRingElem, n_maximal_polyhedra(Sigma))
    TropL = tropical_linear_space(Sigma,multiplicities,convention(nu))
    if !weighted_polyhedral_complex_only
        set_attribute!(TropL,:algebraic_matrix,A)
        set_attribute!(TropL,:tropical_semiring_map,nu)
    end
    return TropL
end

# helper function that constructs a basis of the vanishing set of a linear ideal
# returns a matrix whose rows are the basis vectors
function basis_of_vanishing_of_linear_ideal(I::MPolyIdeal)
    # check that the input ideal is homogeneous and of degree 1
    @req all(isone, total_degree.(Iterators.flatten(monomials.(gens(I))))) "generators of ideal must be linear and homogeneous"

    x = gens(base_ring(I))
    G = gens(I)
    macaulayMatrix = matrix([[coeff(g,xi) for xi in x] for g in G])
    A = transpose(kernel(macaulayMatrix, side = :right))

    return A
end

@doc raw"""
    tropical_linear_space(I::MPolyIdeal, nu::TropicalSemiringMap; weighted_polyhedral_complex_only::Bool=false)

Return the tropicalization of the vanishing set of `I` with respect to the tropical semiring map `nu`.  Requires the generators of `I` to be linear and homogeneous.  If `weighted_polyhedral_complex==true`, will not cache any extra information.

# Examples
```jldoctest
julia> R,(x1,x2,x3,x4) = polynomial_ring(QQ,4);

julia> I = ideal(R,[-x1+x3,-x2+x4])
Ideal generated by
  -x1 + x3
  -x2 + x4

julia> nu = tropical_semiring_map(QQ)
Map into Min tropical semiring encoding the trivial valuation on Rational field

julia> tropical_linear_space(I, nu)
Min tropical linear space

```
"""
function tropical_linear_space(I::MPolyIdeal, nu::Union{Nothing,TropicalSemiringMap}=nothing; weighted_polyhedral_complex_only::Bool=false)
    A = basis_of_vanishing_of_linear_ideal(I)
    TropL = tropical_linear_space(A,nu,
                                  weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
    if !weighted_polyhedral_complex_only
        set_attribute!(TropL,:algebraic_ideal,I)
        set_attribute!(TropL,:tropical_semiring_map,nu)
    end
    return TropL
end



###############################################################################
#
#  Properties
#
###############################################################################
@doc raw"""
    pluecker_indices(TropL::TropicalLinearSpace)

Return the Pluecker indices used to construct `TropL`.
"""
@attr Vector{Vector{Int}} function pluecker_indices(TropL::TropicalLinearSpace)
    # there are two ways to compute the pluecker indices:
    # - from the algebraic matrix
    # - from the tropical matrix
    # (we generally assume that a tropical linear space has at most one of the two)

    # check if tropical matrix is available and, if yes, compute the pluecker indices from it
    if has_attribute(TropL,:tropical_matrix)
        A = tropical_matrix(TropL)
        plueckerIndices, plueckerVector = compute_pluecker_indices_and_vector(A)
        set_attribute!(TropL, :tropical_pluecker_vector, plueckerVector)
        return plueckerIndices
    end

    # otherwise, get the algebraic matrix and compute the pluecker indices from it
    # (will rightfully error if algebraic matrix cannot be found or computed)
    A = algebraic_matrix(TropL)
    plueckerIndices, plueckerVector = compute_pluecker_indices_and_vector(A)
    set_attribute!(TropL, :algebraic_pluecker_vector, plueckerVector)
    return plueckerIndices
end


@doc raw"""
    tropical_pluecker_vector(TropL::TropicalLinearSpace)

Return the tropical Pluecker vector of `TropL`.
"""
@attr Vector function tropical_pluecker_vector(TropL::TropicalLinearSpace)
    A = tropical_matrix(TropL)
    plueckerIndices, plueckerVector = compute_pluecker_indices_and_vector(A)
    set_attribute!(TropL, :pluecker_indices, plueckerIndices)
    return plueckerVector
end


@doc raw"""
    algebraic_pluecker_vector(TropL::TropicalLinearSpace)

Return the Pluecker vector over a valued field used to construct `TropL`.
"""
@attr Vector function algebraic_pluecker_vector(TropL::TropicalLinearSpace)
    A = algebraic_matrix(TropL)
    plueckerIndices, plueckerVector = compute_pluecker_indices_and_vector(A)
    set_attribute!(TropL, :pluecker_indices, plueckerIndices)
    return plueckerVector
end


@doc raw"""
    tropical_semiring_map(TropL::TropicalLinearSpace)

Return the tropical semiring map used to construct `TropL`.  Raises an error, if it is not cached.
"""
function tropical_semiring_map(TropL::TropicalLinearSpace)
    @req has_attribute(TropL,:tropical_semiring_map) "no tropical semiring map cached"
    return get_attribute(TropL, :tropical_semiring_map)::TropicalSemiringMap
end


@doc raw"""
    tropical_matrix(TropL::TropicalLinearSpace)

Return the tropical matrix used to construct `TropL`.  Raises an error, if it is not cached.
"""
function tropical_matrix(TropL::TropicalLinearSpace{minOrMax}) where minOrMax
    @req has_attribute(TropL,:tropical_matrix) "no tropical matrix cached"
    return get_attribute(TropL, :tropical_matrix)::dense_matrix_type(TropicalSemiringElem{minOrMax})
end


@doc raw"""
    algebraic_matrix(TropL::TropicalLinearSpace)

Return the matrix over a valued field used to construct `TropL`.
"""
@attr MatElem function algebraic_matrix(TropL::TropicalLinearSpace)
    I = algebraic_ideal(TropL)
    return basis_of_vanishing_of_linear_ideal(I)
end


@doc raw"""
    algebraic_ideal(TropL::TropicalLinearSpace)

Return the polynomial ideal over a valued field used to construct `TropL`.
Raises an error, if it is not cached.
"""
function algebraic_ideal(TropL::TropicalLinearSpace)
    @req has_attribute(TropL,:algebraic_ideal) "no algebraic ideal cached"
    return get_attribute(TropL, :algebraic_ideal)::MPolyIdeal
end
