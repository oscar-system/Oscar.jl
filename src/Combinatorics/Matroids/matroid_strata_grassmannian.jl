export
    matroid_stratum_matrix_coordinates, matroid_realization_space

@doc Markdown.doc"""
    matroid_stratum_matrix_coordinates(M::Matroid, B::Vector{Int}, F::AbstractAlgebra.Ring = ZZ)

Return the data of the coordinate ring of the matroid stratum of M in the Grassmannian with respect to matrix coordinates. Here, `B` is a basis of `M`` and the submatrix with columns indexed by `B' is the identity. This function returns a pair `(A, W)` where `A` is the coordinate matrix, and `W` is the coordinate ring of the stratum, in general this is a localized quotient ring. 

# Examples
```jldoctest
julia> M = fano_matroid();
julia> (A, W) = matroid_stratum_matrix_coordinates(M, [1,2,4], GF(2));
julia> A # The coordinate matrix with entries in the polynomial ring `R`.
julia> W # The coordinate ring of the stratum; in general a localized quotient ring `(R/I)[S⁻¹]`.
```
"""

function matroid_stratum_matrix_coordinates(M::Matroid, B::GroundsetType,  F::AbstractAlgebra.Ring = ZZ)

    d = rank(M)
    n = length(matroid_groundset(M))

    goodM = isomorphic_matroid(M, [i for i in 1:n])
    #Vector{Int}set difference julia
    goodB = sort!(Int.([M.gs2num[j] for j in B]))

    Bs = bases(goodM)

    goodB in Bs || error("B is not a basis")
    
    R, x, xdict = make_polynomial_ring(Bs,goodB,F)
    return matroid_stratum_matrix_coordinates_given_ring(d, n, goodM, F, goodB, R, x, xdict)
end



@doc Markdown.doc"""
    matroid_realization_space(M::Matroid, A::Vector{Int}, F::AbstractAlgebra.Ring=ZZ)

Returns the data of the coordinate ring of the realization space of
the matroid `M` using matrix coordinates. The matroid `M` should be
a simple and connected matroid, say its rank is ``d``, and ground set
``[n]``. The vector `A` is `rank(M)+1` consists of ``d+1`` elements
(in order) of ``[n]`` such that each ``d``-element subset is a basis of ``M``.

This function returns a pair `(X, W)` where `X` is the
reduced ``d×n`` matrix of variables, and the coordinate ring of
the matroid realization space is `W`.

# Examples
```jldoctest
julia> M = fano_matroid();
julia> (X, W) = matroid_realization_space(M, [1,2,4,7], GF(2));
julia> X # The coordinate matrix.
julia> W # The coordinate ring of the stratum.
```
"""
function matroid_realization_space(M::Matroid, A::GroundsetType, F::AbstractAlgebra.Ring=ZZ)

    n_connected_components(M) == 1 || error("Matroid is not connected")
    
    d = rank(M)
    n = length(matroid_groundset(M))

    goodM = isomorphic_matroid(M, [i for i in 1:n])
    #Vector{Int}
    goodA = sort!(Int.([M.gs2num[j] for j in A]))

    Bs = bases(goodM)

    all([setdiff(goodA,[i]) in Bs for i in goodA]) || error("elements in A are not in general position")
    
    R, x, xdict = realization_polynomial_ring(Bs,goodA,F)
    return matroid_realization_space_given_ring(d, n, goodM, F, goodA, R, x, xdict)
end


# given the bases Bs of a matroid, and a fixed basis B, this function finds
# the nonzero coordinates xij of the coordinate ring of the matroid stratum,
# These correspond to all elements A of Bs such that the symmetric difference
# with B has exactly 2 elements. 

function bases_matrix_coordinates(Bs::Vector{Vector{Int}}, B::Vector{Int})
    
    coord_bases = [b for b in Bs if length(symdiff(B,b)) == 2]
    
    new_coords = Vector{Vector{Int}}([])
    
    for b in coord_bases
        row_b = setdiff(B,b)[1]
        row_b = count(a->(a<row_b), B) + 1
        # count(f, v) does what?
        # - iterate through the elements a in v
        # - compute f(a) 
        # - if that is true, increment the counting variable by 1
        # - otherwise, continue
        # - return the value of the internal counter.
        # Similar with all(a->(a<row_b), B), for instance.
        #row_b = length([a for a in B if a < row_b]) + 1
        
        col_b = setdiff(b,B)[1]
        col_b = col_b - length([a for a in B if a ≤ col_b]) 
        
        push!(new_coords, [row_b,col_b])
    end
    return sort!(new_coords, by = x -> (x[2], x[1]))    
end



# Given the bases Bs of a matroid, a fixed basis B, and a coefficient field F
# this function creates a polynomial ring in xij, where the xij are 
# determined by the function basis_matrix_coordinates. This function also
# returns (as the 3rd element of a triple) a dictionary (i,j) => xij.

function make_polynomial_ring(Bs::Vector{Vector{Int}}, B::Vector{Int},
                              F::AbstractAlgebra.Ring)
    
    MC = bases_matrix_coordinates(Bs, B)
    R, x = PolynomialRing(F, :"x"=>MC)
    xdict = Dict{Vector{Int}, MPolyElem}([MC[i] => x[i] for i in 1:length(MC)])
    return R, x, xdict
end


# This function returns a d x (n-d) matrix with values in the polynomial ring
# created by make_polynomial_ring. The entries are xij, except where the
# value is 0, as determined by the nonbases. 

function make_coordinate_matrix_no_identity(d::Int, n::Int,
                                            MC::Vector{Vector{Int}},
                                            R::MPolyRing, x::Vector{T},
                                            xdict::Dict{Vector{Int}, MPolyElem}) where T <: MPolyElem
    
    S = MatrixSpace(R, d, n-d)
    X = S()
    for j in 1:n-d, i in 1:d
        if [i,j] in MC
            X[i,j] = xdict[[i,j]]
        else
            X[i,j] = R(0)
        end
    end
    return X
end


# M and N have same number of rows, and M has #B columns, both have entries in ring R
function interlace_columns(M, N, B::Vector{Int},
                           R::MPolyRing, x::Vector{T}) where T <: MPolyElem 
    
    M_nrows, M_ncols = size(M)
    N_nrows, N_ncols = size(N)
    n = M_ncols + N_ncols

    S = MatrixSpace(R, M_nrows, n)
    Bc = [i for i in 1:n if !(i in B)]
    
    X = S()
    X[:, B] = M
    X[:, Bc] = N
    
    return X 
end

# This makes the matrix X from which we compute the coordinate ring of the matroid
# stratum. It has the identity matrix at columns indexed by B, 0's at locations
# determined by the nonbases of X. 
function make_coordinate_matrix(d::Int, n::Int, MC::Vector{Vector{Int}},
                                B::Vector{Int},
                                R::MPolyRing, x::Vector{T},
                                xdict::Dict{Vector{Int}, MPolyElem}) where T <: MPolyElem
    
    Id = identity_matrix(R,d)
    Xpre = make_coordinate_matrix_no_identity(d, n, MC, R, x, xdict)
    return interlace_columns(Id, Xpre, B, R, x)
end


# This function returns all d x d determinants of the matrix X from above
# of all collections of d-columns coming from the bases of the matroid.

function bases_determinants(X::MatrixElem{T}, Bs::Vector{Vector{Int}}) where {T<:MPolyElem}
    return unique!([det(X[:, b]) for b in Bs ])
end


#function bases_determinants(X::Matrix{T}, Bs::Vector{Vector{Int}})  where {T<:MPolyElem}

    #d::Int, n::Int, Bs::Vector{Vector{Int}},
    #MC::Vector{Vector{Int}},
    #B::Vector{Int}, R::MPolyRing, x::Vector{T},
    #xdict::Dict{Vector{Int}, MPolyElem}) where T <: MPolyElem
    
    #X = make_coordinate_matrix(d, n, MC, B, R, x, xdict)
    
#    return unique!([det(X[:, b]) for b in Bs ])
#end


# This forms the semigroup of the polynomial ring from make_polynomial_ring
# generated by the determinants from bases_determinants.
# function localizing_semigroup(d::Int, n::Int, Bs::Vector{Vector{Int}},
#                               MC::Vector{Vector{Int}}, B::Vector{Int},
#                               R::MPolyRing, x::Vector{T},
#                               xdict::Dict{Vector{Int}, MPolyElem}) where T <: MPolyElem

    
#     basesX = bases_determinants(X,Bs)
#     #basesX = bases_determinants(d, n, Bs, MC, B, R, x, xdict)
    
#     sTotal = MPolyPowersOfElement(basesX[1])
    
#     if length(basesX) == 1
#         return sTotal
#     end
    
#     for i in 2:length(basesX)
#         if !(basesX[i] in sTotal)
#             sTotal = product(sTotal, MPolyPowersOfElement(basesX[i]))
#         end
#     end
    
#     return sTotal
    
# end


# This function returns a triple (X, SinvR, Ipre) where X is the matrix
# from make_coordinate_matrix, SinvR is the localization of the
# polynomial ring at the semigroup created by localizing_semigroup,
# and Ipre is the ideal generated by the d x d determinants of nonbases.
# Note that Ipre is not saturated. 
function matroid_stratum_matrix_coordinates_given_ring(d::Int, n::Int,
                                                       M::Matroid,
                                                       F::AbstractAlgebra.Ring,
                                                       B::Vector{Int},
                                                       R::MPolyRing,
                                                       x::Vector{T},
                                                       xdict::Dict{Vector{Int}, MPolyElem}) where T <: MPolyElem
    

    Bs = bases(M)
    NBs = nonbases(M)
    NBsNotVariable = [nb for nb in NBs if length(symdiff(B,nb)) != 2]

    MC = bases_matrix_coordinates(Bs,B)

    
    X =  make_coordinate_matrix(d, n, MC, B, R, x, xdict)                        
    basesX = bases_determinants(X, Bs)
    
    #S = localizing_semigroup(d, n, Bs, MC, B, R, x, xdict)
    
    S = MPolyPowersOfElement(R , basesX)
    SinvR , iota = Localization(R, S)
    # X = make_coordinate_matrix(d, n, MC, B, R, x, xdict)

    Igens = unique!([det(X[:, nb]) for nb in NBsNotVariable ])
    Iloc = ideal(SinvR, Igens)
    if iszero(Iloc)
      return (X, SinvR)
    else
      W, _ = quo(SinvR, Iloc)
      return (X, W)
    end

end

    
function realization_bases_coordinates(Bs::Vector{Vector{Int}}, A::Vector{Int})

    d = length(Bs[1])
    B = A[1:d]
    c1 = A[d+1]

    coord_bases = [b for b in Bs if length(symdiff(B,b)) == 2]
    
    new_coords = Vector{Vector{Int}}()
    
    for b in coord_bases
        if is_subset(b, A)
            continue
        end
        
        row_b = setdiff(B,b)[1]
        row_b = count(a->(a<row_b), B) + 1
        
        col_b = setdiff(b,B)[1]
        col_b = col_b - length([a for a in A if a <= col_b])

        push!(new_coords, [row_b, col_b])
    end
    
    return sort!(new_coords, by = x -> (x[2], x[1]))
end


function partial_matrix_max_rows(Vs::Vector{Vector{Int}})
    nr = maximum([x[1] for x in Vs])
    cols = unique!([x[2] for x in Vs])
    first_nonzero_cols = Dict{Int, Int}(c => maximum(i for i in 1:nr if [i,c] in Vs) for c in cols)
    return first_nonzero_cols 
end



function realization_polynomial_ring(Bs::Vector{Vector{Int}}, A::Vector{Int},
                                     F::AbstractAlgebra.Ring)
    
    MC = realization_bases_coordinates(Bs, A)
    D = partial_matrix_max_rows(MC)
    MR = [x for x in MC if x[1] != D[x[2]]]
    R, x = PolynomialRing(F, :"x"=>MR)
    xdict = Dict{Vector{Int}, MPolyElem}(MR[i] => x[i] for i in 1:length(MR))
    return R, x, xdict
end


function matrix_realization_small(d::Int, n::Int, MC::Vector{Vector{Int}},
                                  R::MPolyRing, x::Vector{T},
                                  xdict::Dict{Vector{Int}, MPolyElem}) where T <: MPolyElem

    D = partial_matrix_max_rows(MC)
    MR = [x for x in MC if x[1] != D[x[2]]]
    S = MatrixSpace(R, d, n-d-1)
    X = S()
    for j in 1:n-d-1, i in 1:d
        if [i,j] in MR
            X[i,j] = xdict[[i,j]]
        elseif(j in keys(D) && i == D[j])
            X[i,j] = R(1)            
        else
            X[i,j] = R(0)
        end
    end
    return X
end

function projective_identity(d::Int)
    if d == 1
        return [1]
    end 
    X = zeros(Int, d, d+1)
    for i in 1:d
        X[i,i] = 1
        X[i,d+1] = 1
    end
    return X
end

function realization_coordinate_matrix(d::Int, n::Int, MC::Vector{Vector{Int}},
                                       A::Vector{Int}, R::MPolyRing, x::Vector{T},
                                       xdict::Dict{Vector{Int}, MPolyElem}) where T <: MPolyElem
    
    Id = projective_identity(d) 
    Xpre = matrix_realization_small(d, n, MC, R, x, xdict)
    return interlace_columns(Id, Xpre, A, R, x)
end



function realization_bases_determinants(X::MatrixElem{T}, Bs::Vector{Vector{Int}}) where {T<:MPolyElem}
    return unique!([det(X[:, b]) for b in Bs ])
end

# For a given list `l` of polynomials in a ring `R` this generates 
# the semigroup generated by all products of powers of elements in `l`.
#
# TODO: This function seems deprecated.
function realization_localizing_semigroup(basesX::Vector{T}) where {T<:MPolyElem}
  iszero(length(basesX)) && error("list must not be empty")
  R = parent(first(basesX))
  all(x->(parent(x)===R), basesX) || error("elements must belong to the same ring")

  return MPolyPowersOfElement(R, basesX)
end


function matroid_realization_space_given_ring(d::Int, n::Int, M::Matroid,
                                              F::AbstractAlgebra.Ring, A::Vector{Int},
                                              R::MPolyRing, x::Vector{T},
                                              xdict::Dict{Vector{Int}, MPolyElem}) where T <: MPolyElem
    
    Bs = bases(M)
    NBs = nonbases(M)
    MC = realization_bases_coordinates(Bs,A)
    NBsNotVariable = [nb for nb in NBs if length(symdiff(A[1:d],nb)) != 2]


    X = realization_coordinate_matrix(d, n, MC, A, R, x, xdict)
    basesX = realization_bases_determinants(X, Bs)
    
    S = MPolyPowersOfElement(R , basesX)
    #S = realization_localizing_semigroup(basesX); 
    SinvR , iota = Localization(R, S)
    

    Igens = [det(X[:, nb]) for nb in NBsNotVariable ]
    Iloc = ideal(SinvR, Igens)
    if iszero(Iloc)
      return (X, SinvR)
    else
      W, _ = quo(SinvR, Iloc)
      return (X, W)
    end
end


