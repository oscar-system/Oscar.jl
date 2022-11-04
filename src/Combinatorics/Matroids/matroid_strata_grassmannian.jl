export
    matroid_stratum_matrix_coordinates


function bases_matrix_coordinates(Bs::Vector{Vector{Int64}}, B::Vector{Int64})
    
    coordBases = [b for b in Bs if length(symdiff(B,b)) == 2]
    
    newCoords = Vector{Vector{Int64}}([])
    
    for b in coordBases
        row_b = setdiff(B,b)[1]
        row_b = length([a for a in B if a < row_b]) + 1
        
        col_b = setdiff(b,B)[1]
        col_b = col_b - length([a for a in B if a ≤ col_b]) 
        
        push!(newCoords, [row_b,col_b])
    end
    return sort!(newCoords, by = x -> (x[2], x[1]))    
end


function make_polynomial_ring(Bs::Vector{Vector{Int64}}, B::Vector{Int64},
                              F::AbstractAlgebra.Ring)
    
    MC = bases_matrix_coordinates(Bs, B)
    R, x = PolynomialRing(F, :"x"=>MC)
    W = ones(Int64, length(MC))
    R, x = grade(R,W)
    xdict = Dict{Vector{Int64}, MPolyElem}([MC[i] => x[i] for i in 1:length(MC)])
    return R, x, xdict
end

function make_coordinate_matrix_no_identity(d::Int64, n::Int64,
                                            MC::Vector{Vector{Int64}},
                                            R::MPolyRing, x::Vector{T},
                                            xdict::Dict{Vector{Int64}, MPolyElem}) where T <: MPolyElem
    
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
function interlace_columns(M, N, B::Vector{Int64},
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


function make_coordinate_matrix(d::Int64, n::Int64, MC::Vector{Vector{Int64}},
                                B::Vector{Int64},
                                R::MPolyRing, x::Vector{T},
                                xdict::Dict{Vector{Int64}, MPolyElem}) where T <: MPolyElem
    
    Id = identity_matrix(R,d)
    Xpre = make_coordinate_matrix_no_identity(d, n, MC, R, x, xdict)
    return interlace_columns(Id, Xpre, B, R, x)
end



function bases_determinants(d::Int64, n::Int64, Bs::Vector{Vector{Int64}},
                            MC::Vector{Vector{Int64}},
                            B::Vector{Int64}, R::MPolyRing, x::Vector{T},
                            xdict::Dict{Vector{Int64}, MPolyElem}) where T <: MPolyElem
    
    X = make_coordinate_matrix(d, n, MC, B, R, x, xdict)
    return unique!([det(X[:, b]) for b in Bs ])
end



function localizing_semigroup(d::Int64, n::Int64, Bs::Vector{Vector{Int64}},
                              MC::Vector{Vector{Int64}}, B::Vector{Int64},
                              R::MPolyRing, x::Vector{T},
                              xdict::Dict{Vector{Int64}, MPolyElem}) where T <: MPolyElem
    
    basesX = bases_determinants(d, n, Bs, MC, B, R, x, xdict)
    
    sTotal = MPolyPowersOfElement(basesX[1])
    
    if length(basesX) == 1
        return sTotal
    end
    
    for i in 2:length(basesX)
        if !(basesX[i] in sTotal)
            sTotal = product(sTotal, MPolyPowersOfElement(basesX[i]))
        end
    end
    
    return sTotal
    
end


function matroid_stratum_matrix_coordinates_given_ring(d::Int64, n::Int64,
                                                       M::Matroid,
                                                       F::AbstractAlgebra.Ring,
                                                       B::Vector{Int64},
                                                       R::MPolyRing,
                                                       x::Vector{T},
                                                       xdict::Dict{Vector{Int64}, MPolyElem}) where T <: MPolyElem
    
    Bs = bases(M)
    NBs = nonbases(M)
    MC = bases_matrix_coordinates(Bs,B)


    S = localizing_semigroup(d, n, Bs, MC, B, R, x, xdict)
    SinvR , iota = Localization(R, S)
    X = make_coordinate_matrix(d, n, MC, B, R, x, xdict)

    Igens = unique!([det(X[:, nb]) for nb in NBs ])    
    Ipre = ideal(SinvR, Igens)
    I = saturated_ideal(Ipre)
    
    return (SinvR,I)
end

@doc Markdown.doc"""
    matroid_stratum_matrix_coordinates(M::Matroid, B::Vector{Int64}, F::AbstractAlgebra.Ring = ZZ)

Return the data of the coordinate ring of the matroid stratum of M in the Grassmannian with respect to matrix coordinates. Here, ```B``` is a basis of ```M``` and the submatrix with columns indexed by ```B''' is the identity. This function returns a pair ```(SinvR, I)``` where the coordinate ring is isomorphic to ```SinvR / I```.

# Examples
```jldoctest
julia> M = fano_matroid()
julia> matroid_stratum_matrix_coordinates(M, [1,2,4], GF(2))
(localization of Multivariate Polynomial Ring in 9 variables x[2, 1], x[1, 1], x[1, 2], x[3, 2], ..., x[3, 4] over Integer Ring graded by 
  x[2, 1] -> [1]
  x[1, 1] -> [1]
  x[1, 2] -> [1]
  x[3, 2] -> [1]
  x[2, 3] -> [1]
  x[3, 3] -> [1]
  x[2, 4] -> [1]
  x[1, 4] -> [1]
  x[3, 4] -> [1] at the powers of MPolyElem_dec{fmpz, fmpz_mpoly}[x[3, 3], x[2, 1]*x[3, 3]*x[1, 4] + x[1, 1]*x[2, 3]*x[3, 4] - x[1, 1]*x[3, 3]*x[2, 4], -x[2, 3]*x[1, 4], x[1, 2]*x[2, 3]*x[3, 4] - x[1, 2]*x[3, 3]*x[2, 4] - x[3, 2]*x[2, 3]*x[1, 4], -x[3, 2], -x[2, 1]*x[1, 2]*x[3, 4] + x[2, 1]*x[3, 2]*x[1, 4] - x[1, 1]*x[3, 2]*x[2, 4], x[1, 2]*x[2, 4], x[2, 1], -x[1, 1]*x[3, 4]], ideal(2, x[1, 2]*x[3, 4] + x[3, 2]*x[1, 4], -x[2, 3]*x[3, 4] + x[3, 3]*x[2, 4], -x[2, 1]*x[1, 4] + x[1, 1]*x[2, 4], x[2, 1]*x[3, 3]*x[1, 4] - x[1, 1]*x[2, 3]*x[3, 4], x[2, 1]*x[1, 2]*x[3, 3] + x[1, 1]*x[3, 2]*x[2, 3]))
```
"""

function matroid_stratum_matrix_coordinates(M::Matroid, B::Vector{Int64},
                                            F::AbstractAlgebra.Ring = ZZ)
    d = rank(M)
    n = length(matroid_groundset(M))
    Bs = bases(M)
    R, x, xdict = make_polynomial_ring(Bs,B,F)
    return matroid_stratum_matrix_coordinates_given_ring(d, n, M, F, B, R, x, xdict)
end
