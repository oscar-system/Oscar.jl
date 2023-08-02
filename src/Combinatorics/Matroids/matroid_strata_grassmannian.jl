struct MatroidRealizationSpace
    defining_ideal::Union{Ideal, NumFieldOrdIdl, Nothing}
    inequations::Union{Vector{Oscar.RingElem},Nothing}
    ambient_ring::Union{Oscar.MPolyRing, Ring, Nothing}
    representation_matrix::Union{Oscar.MatElem,Nothing}
    representable::Bool
    F::AbstractAlgebra.Ring
    char::Union{Int,Nothing}
    q::Union{Int,Nothing}
end

function Base.show(io::IO, RS::MatroidRealizationSpace)
    if !RS.representable
        print(io, "The matroid is not representable.")
    else
        println(io, "The representations are parametrized by")
        # println isn't ideal as it prints the matrix as one big line
        display(RS.representation_matrix)
        println(io, "in the ", RS.ambient_ring)
        I = RS.defining_ideal
        if (typeof(I) <: NumFieldOrdIdl && I.gen != ZZ(0)) || (typeof(I) <: Ideal && !iszero(I))
            println(io, "within the vanishing set of the ideal\n", RS.defining_ideal)
        end
        if length(RS.inequations) > 0
            println(io, "avoiding the zero loci of the polynomials\n", RS.inequations)
        end
    end
end

function is_representable(M::Matroid; char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing)::Bool
    RS = realization_space(M, char=char, q=q)
    return RS.representable
end

function defining_ideal(RS::MatroidRealizationSpace)
    return RS.defining_ideal
end

function inequations(RS::MatroidRealizationSpace)
    return RS.inequations
end

function ambient_ring(RS::MatroidRealizationSpace)
    return RS.ambient_ring
end

function representation_matrix(RS::MatroidRealizationSpace)
    return RS.representation_matrix
end

#=

@doc raw"""
    matroid_stratum_matrix_coordinates(M::Matroid, B::GroundsetType, F::AbstractAlgebra.Ring = ZZ)

Return the data of the coordinate ring of the matroid stratum of M in the Grassmannian with respect to matrix coordinates. Here, `B` is a basis of `M`` and the submatrix with columns indexed by `B' is the identity. This function returns a pair `(A, W)` where `A` is the coordinate matrix, and `W` is the coordinate ring of the stratum, in general this is a localized quotient ring. 

# Examples
```jldoctest
julia> M = fano_matroid()

julia> (A, W) = matroid_stratum_matrix_coordinates(M, [1,2,4], GF(2))

julia> A # The coordinate matrix with entries in the polynomial ring `R`.
[1   0   x[1, 1]   0   x[1, 2]         0   x[1, 4]]
[0   1   x[2, 1]   0         0   x[2, 3]   x[2, 4]]
[0   0         0   1   x[3, 2]   x[3, 3]   x[3, 4]]

julia> W # The coordinate ring of the stratum in general a localized quotient ring `(R/I)[S⁻¹]`.
Localization of Quotient of Multivariate Polynomial Ring in 9 variables x[1, 1], x[2, 1], x[1, 2], x[3, 2], ..., x[3, 4] over Galois field with characteristic 2 by ideal(x[2, 3]*x[3, 4] + x[3, 3]*x[2, 4], x[1, 2]*x[3, 4] + x[3, 2]*x[1, 4], x[1, 1]*x[2, 4] + x[2, 1]*x[1, 4], x[1, 1]*x[3, 2]*x[2, 3] + x[2, 1]*x[1, 2]*x[3, 3]) at the multiplicative set powers of fpRingElem[x[3, 3]*x[1, 4], x[1, 1]*x[2, 3]*x[3, 4] + x[1, 1]*x[3, 3]*x[2, 4] + x[2, 1]*x[3, 3]*x[1, 4], x[2, 3]*x[1, 4], x[1, 2]*x[2, 3]*x[3, 4] + x[1, 2]*x[3, 3]*x[2, 4] + x[3, 2]*x[2, 3]*x[1, 4], x[3, 2]*x[2, 4], x[1, 1]*x[3, 2]*x[2, 4] + x[2, 1]*x[1, 2]*x[3, 4] + x[2, 1]*x[3, 2]*x[1, 4], x[1, 2]*x[2, 4], x[2, 4], x[1, 4], x[2, 1]*x[3, 4], x[1, 1]*x[3, 4], x[3, 4], x[3, 2]*x[2, 3], x[1, 2]*x[3, 3], x[1, 2]*x[2, 3], x[2, 3], x[1, 1]*x[2, 3], x[2, 1]*x[3, 3], x[1, 1]*x[3, 3], x[3, 3], x[1, 2], x[2, 1]*x[1, 2], x[2, 1]*x[3, 2], x[1, 1]*x[3, 2], x[3, 2], x[2, 1], x[1, 1], 1]
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



@doc raw"""
    matroid_realization_space(M::Matroid, A::GroundsetType, F::AbstractAlgebra.Ring=ZZ)

Return the data of the coordinate ring of the realization space of
the matroid `M` using matrix coordinates. The matroid `M` should be
a simple and connected matroid, say its rank is ``d``, and ground set
``[n]``. The vector `A` is `rank(M)+1` consists of ``d+1`` elements
(in order) of ``[n]`` such that each ``d``-element subset is a basis of ``M``.

This function returns a pair `(X, W)` where `X` is the
reduced ``d×n`` matrix of variables, and the coordinate ring of
the matroid realization space is `W`.

# Examples
```jldoctest
julia> M = fano_matroid()

julia> (X, W) = matroid_realization_space(M, [1,2,4,7], GF(2))

julia> X # The coordinate matrix.
[1   0   x[1, 1]   0   x[1, 2]         0   1]
[0   1         1   0         0   x[2, 3]   1]
[0   0         0   1         1         1   1]

julia> W # The coordinate ring of the stratum.
Localization of Quotient of Multivariate Polynomial Ring in x[1, 1], x[1, 2], x[2, 3] over Galois field with characteristic 2 by ideal(x[2, 3] + 1, x[1, 2] + 1, x[1, 1] + 1, x[1, 1]*x[2, 3] + x[1, 2]) at the multiplicative set powers of fpRingElem[1, x[1, 1]*x[2, 3] + x[1, 1] + 1, x[2, 3], x[1, 2]*x[2, 3] + x[1, 2] + x[2, 3], x[1, 1] + x[1, 2] + 1, x[1, 2], x[1, 1], x[1, 2]*x[2, 3], x[1, 1]*x[2, 3]]
```
"""
function matroid_realization_space(M::Matroid, A::GroundsetType, F::AbstractAlgebra.Ring=ZZ)

    n_connected_components(M) == 1 || error("Matroid is not connected")
    is_simple(M) || error("Matroid is not simple")
    
    d = rank(M)
    n = length(matroid_groundset(M))

    if d == 1
        return F
    end

    
    goodM = isomorphic_matroid(M, [i for i in 1:n])
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
    R, x = polynomial_ring(F, :"x"=>MC)
    xdict = Dict{Vector{Int}, RingElem}([MC[i] => x[i] for i in 1:length(MC)])
    return R, x, xdict
end


# This function returns a d x (n-d) matrix with values in the polynomial ring
# created by make_polynomial_ring. The entries are xij, except where the
# value is 0, as determined by the nonbases. 

function make_coordinate_matrix_no_identity(d::Int, n::Int,
                                            MC::Vector{Vector{Int}},
                                            R::MPolyRing, x::Vector{T},
                                            xdict::Dict{Vector{Int}, RingElem}) where T <: RingElem
    
    X = zero_matrix(R, d, n-d)
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
function interlace_columns(M::MatrixElem{T}, N::MatrixElem{T}, B::Vector{Int},
                           R::MPolyRing, x::Vector{T}) where T <: RingElem 
    
    M_nrows, M_ncols = size(M)
    N_nrows, N_ncols = size(N)
    n = M_ncols + N_ncols

    Bc = [i for i in 1:n if !(i in B)]
    
    X = zero_matrix(R, M_nrows, n)
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
                                xdict::Dict{Vector{Int}, RingElem}) where T <: RingElem
    
    Id = identity_matrix(R,d)
    Xpre = make_coordinate_matrix_no_identity(d, n, MC, R, x, xdict)
    return interlace_columns(Id, Xpre, B, R, x)
end


# This function returns all d x d determinants of the matrix X from above
# of all collections of d-columns coming from the bases of the matroid.

function bases_determinants(X::MatrixElem{T}, Bs::Vector{Vector{Int}}) where {T<:RingElem}
    return unique!([det(X[:, b]) for b in Bs ])
end


#function bases_determinants(X::Matrix{T}, Bs::Vector{Vector{Int}})  where {T<:RingElem}

    #d::Int, n::Int, Bs::Vector{Vector{Int}},
    #MC::Vector{Vector{Int}},
    #B::Vector{Int}, R::MPolyRing, x::Vector{T},
    #xdict::Dict{Vector{Int}, RingElem}) where T <: RingElem
    
    #X = make_coordinate_matrix(d, n, MC, B, R, x, xdict)
    
#    return unique!([det(X[:, b]) for b in Bs ])
#end


# This forms the semigroup of the polynomial ring from make_polynomial_ring
# generated by the determinants from bases_determinants.
# function localizing_semigroup(d::Int, n::Int, Bs::Vector{Vector{Int}},
#                               MC::Vector{Vector{Int}}, B::Vector{Int},
#                               R::MPolyRing, x::Vector{T},
#                               xdict::Dict{Vector{Int}, RingElem}) where T <: RingElem

    
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


# This function returns the output that appears in matroid_stratum_matrix_coordinates. 
function matroid_stratum_matrix_coordinates_given_ring(d::Int, n::Int,
                                                       M::Matroid,
                                                       F::AbstractAlgebra.Ring,
                                                       B::Vector{Int},
                                                       R::MPolyRing,
                                                       x::Vector{T},
                                                       xdict::Dict{Vector{Int}, RingElem}) where T <: RingElem
    

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
    R, x = polynomial_ring(F, :"x"=>MR)
    xdict = Dict{Vector{Int}, RingElem}(MR[i] => x[i] for i in 1:length(MR))
    return R, x, xdict
end


function matrix_realization_small(d::Int, n::Int, MC::Vector{Vector{Int}},
                                  R::MPolyRing, x::Vector{T},
                                  xdict::Dict{Vector{Int}, RingElem}) where T <: RingElem

    D = partial_matrix_max_rows(MC)
    MR = [x for x in MC if x[1] != D[x[2]]]
    X = zero_matrix(R, d, n-d-1)
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
        return ones(Int, 1, 1)
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
                                       xdict::Dict{Vector{Int}, RingElem}) where T <: RingElem
    
    Id = matrix(R, projective_identity(d))
    Xpre = matrix_realization_small(d, n, MC, R, x, xdict)
    return interlace_columns(Id, Xpre, A, R, x)

end



function realization_bases_determinants(X::MatrixElem{T}, Bs::Vector{Vector{Int}}) where {T<:RingElem}
    return unique!([det(X[:, b]) for b in Bs ])
end




function matroid_realization_space_given_ring(d::Int, n::Int, M::Matroid,
                                              F::AbstractAlgebra.Ring, A::Vector{Int},
                                              R::MPolyRing, x::Vector{T},
                                              xdict::Dict{Vector{Int}, RingElem}) where T <: RingElem
    
    Bs = bases(M)
    NBs = nonbases(M)
    MC = realization_bases_coordinates(Bs,A)
    NBsNotVariable = [nb for nb in NBs if length(symdiff(A[1:d],nb)) != 2]


    X = realization_coordinate_matrix(d, n, MC, A, R, x, xdict)
    basesX = realization_bases_determinants(X, Bs)
    
    S = MPolyPowersOfElement(R , basesX)
    #S = realization_localizing_semigroup(basesX) 
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
=#



function realization_space_matrix(M::Matroid, B::Vector{Int}, F::AbstractAlgebra.Ring)
    # prepare the combinatorial data
    
    circs = fundamental_circuits_of_basis(M,B)
    
    nonIdCols = setdiff(matroid_groundset(M),B)
    circs = [setdiff( c, nonIdCols ) for c in circs]

    rk = rank(M)
    n = length(M)

    # we start by computing the number of variables:
    numVars = sum([length(intersect(c,B))-1 for c in circs])-(rk-1)+(length(connected_components(M))-1)
    
    if numVars > 0
        R, x = polynomial_ring(F, numVars)
    else
        R = F
        x = Vector{MPolyRingElem}()
    end

    unUsedRowsForOnes = collect(2:rk)
    
    # create the matrix and fill it with entries
    
    mat = zero_matrix(R,rk,n)
    
    for i in 1:rk
        mat[i,B[i]]=R(1)
    end   
    
    varCounter = 1
    
    for col in 1:(n-rk)
        for row in 1:rk
            circ = circs[col]
            c = nonIdCols[col]
            
            if B[row] == minimum(circ)
                mat[row,c] = R(1)
            elseif B[row] in circ 
                if row in unUsedRowsForOnes
                    mat[row,c] = R(1)
                    unUsedRowsForOnes = setdiff(unUsedRowsForOnes,[row])
                else
                    mat[row,c] = x[varCounter]
                    varCounter = varCounter+1
                end
            else
                mat[row, c] = R(0)
            end
        end
    end
    
    return (R, x, mat)
end


# Given a matroid M with a basis B, this functions computes for a i in groundset(M)\B a circuit in B \cup i.
# This is a function used in the construction of the matrix underlying the representation space computation. 
function fundamental_circuits_of_basis(M::Matroid, B::Vector{Int})
    remaining_elts = setdiff(matroid_groundset(M),B)
    fund_circs = Vector{Vector{Int}}()
    for i in remaining_elts
         push!(fund_circs,fundamental_circuit(M,B,i))
     end
    return fund_circs
 end

 function realization_space(M::Matroid; B::Union{GroundsetType,Nothing} = nothing, 
    F::AbstractAlgebra.Ring = ZZ, saturate::Bool=false, 
    char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing)::MatroidRealizationSpace

    if char!=nothing && !isprime(char) && char!=0
        error("The characteristic has to be 0 or a prime number.")
    end

    #Construct the base ring as F_p if q=p^k
    if q!=nothing
        isprimepower, p, k = is_prime_power_with_data(q)
        if !isprimepower
            error("The given q has to be a prime power.")
        end
        if char!=nothing && char!=p
            error("The given characteristic doesn't match q.")
        else
            char = p
        end
    end

    if char == 0
        F = QQ
    elseif char != nothing
        F = GF(char)
    end

    rk = rank(M)
    n = length(M)

    goodM = isomorphic_matroid(M, [i for i in 1:n])

    Bs = bases(goodM)
    
    if !isnothing(B) 
        goodB = sort!(Int.([M.gs2num[j] for j in B]))
    else
        goodB = find_good_basis_heuristically(goodM)
    end

    polyR, x, mat = realization_space_matrix(goodM, goodB, F)
    
    eqs = Vector{RingElem}()
    ineqs = Vector{RingElem}()

    #need to catch corner-case if there are no variables at all
    if length(x) == 0
        return MatroidRealizationSpace(ideal(polyR,0), ineqs, polyR, mat, true, F, char, q)
    end
    
    for col in subsets(Vector(1:n),rk)
        
        col_det = det(mat[:,col])
        
        if total_degree( col_det ) <= 0 
            
            if col_det != 0 && col in Bs 
                if isunit( col_det ) 
                    continue
                end
            elseif col_det != 0 # and col is not a basis
                erorr("determinant nonzero but set not a basis")
            elseif col in Bs 
                error( "determinant zero but set is a basis" )
            else
                continue
            end
            
        end
        
        if  col in Bs
            push!( ineqs, col_det )
        else
            push!( eqs, col_det )
        end
        
    end

    if q != nothing
        for x_elem in x
            push!(eqs, x_elem^q-x_elem)
        end
    end

    def_ideal = ideal(polyR,eqs)
    def_ideal = ideal(groebner_basis(def_ideal))
#    if !iszero(def_ideal)
#        for i in 1:length(ineqs)
#            ineqs[i] = reduce(ineqs[i], gens(def_ideal))
#        end
#    end
    ineqs = gens_2_factors(ineqs)
    
    if saturate || polyR.nvars < 10
        def_ideal = stepwise_saturation(def_ideal,ineqs)
        def_ideal = ideal(groebner_basis(def_ideal))
    end
        
    representable = !(isone(def_ideal))

    !representable && return MatroidRealizationSpace(def_ideal, nothing, nothing, nothing, false, F, char, q)

#    if !representable
#        return MatroidRealizationSpace(def_ideal, nothing, nothing, nothing, false)
#    end

    return MatroidRealizationSpace(def_ideal, ineqs, polyR, mat, representable, F, char, q)
end


# A heuristic function that tries to find a sensible basis for the moduli space computation for which the defining ideal is not too complicated
function find_good_basis_heuristically(M::Matroid)
    bs = bases(M)
    cs = circuits(M)
    min_num_vars = length(cs)*rank(M)
    min_basis = bs[1]
    for bi in 1:length(bs)
        current_num_vars = 0
        for c in cs
            for e in c
                @inbounds if e in bs[bi]
                    current_num_vars += 1
                end
            end
        end
        if current_num_vars < min_num_vars
            min_num_vars = current_num_vars
            min_basis = bs[bi]
        end
    end
    return min_basis
end


#returns the factors of f, but not the exponents
function poly_2_factors(f::RingElem)
    return collect(keys(Dict(factor(f))))
end

# returns the unique factors of the elements of Sgen, again no exponents. 
function gens_2_factors(Sgens::Vector{<:RingElem})
    return unique!(vcat([poly_2_factors(f) for f in Sgens]...))
end


function stepwise_saturation(I::MPolyIdeal, Sgens::Vector{<:RingElem})
    foreach(f -> I = saturation(I,ideal([f])), Sgens)
    return I
end








#####################
# initial reduction #
#####################


function update_hom(v_elim::RingElem, w_repl::RingElem, phi::Oscar.MPolyAnyMap)
    R = codomain(phi)
    phi2_im = replace(gens(R), v_elim => w_repl)
    phi2 = hom(R, R, a->a, phi2_im)
    return hom(R, R, a->a, phi2.(phi.(gens(R))))
end

function small_find_solution(v::RingElem, Igens::Vector{<:RingElem})
    
#    v_deg_1_no_coef = [g for g in Igens if isone(degree(g,v)) && is_unit(coeff(g,[v],[1])) ]
#    println("v_deg_1_no_coef: ", v_deg_1_no_coef)
#    length(v_deg_1_no_coef) != 0 || return "can't isolate"
#    f = first(v_deg_1_no_coef)
#    return -coeff(f, [v], [0]) / coeff(f, [v], [1])
        
    with_v_deg_1 = [g for g in Igens if (total_degree(g) == 1 && is_unit(coeff(g,[v],[1])))] 
    length(with_v_deg_1) == 0 && return "can't isolate"
    
    f = first(with_v_deg_1)
    return -coeff(f, [v], [0])/coeff(f,[v],[1])
end


function small_reduce_one_step(Igens::Vector{<:RingElem}, phi::Oscar.MPolyAnyMap, 
                               elim::Vector{<:RingElem}, fullyReduced::Bool)
    Ivars = ideal_vars(Igens); 
    length(Igens) == 0 && return (Igens, phi, elim, true)
    for v in Ivars 
        w = small_find_solution(v, Igens)
        w isa String && continue
        
        push!(elim, v)
        phi = update_hom(v, w, phi)
        Igens = filter(x->x!=0, phi.(Igens))
        return (Igens, phi, elim, false)
    end
    return (Igens, phi, elim, true)
end


function small_reduce_full_rec(Igens::Vector{<:RingElem}, 
                               phi::Union{<:Oscar.MPolyAnyMap,Nothing} = nothing, 
                               elim::Vector{<:RingElem} = Vector{RingElem}(), 
                               fullyReduced=false)

    if isnothing(phi)
        R = parent(Igens[1])
        phi=hom(R, R, a->a, gens(R))
    end

    output = small_reduce_one_step(Igens, phi, elim, fullyReduced)
    output isa String && return "Not Realizable 0 in Semigroup"
    
    (Igens, phi, elim, fullyReduced) = output    
    !fullyReduced && return small_reduce_full_rec(Igens, phi, elim, fullyReduced)
    
    R = domain(phi)
    x = gens(R)    
    cR = coefficient_ring(R)
    
    if length(elim) == length(x)
        S = cR
    else
        S, y = polynomial_ring(cR, length(x) - length(elim))
    end        
    z=[]
    j=1
    for ele in x
        if ele in elim
            push!(z,0)
        else
            push!(z, y[j])
            j+=1
        end
    end
    phi2 = hom(R, S, a->a, z)
    phif = phi*phi2

    return phif
end


function small_reduce(MRS::MatroidRealizationSpace) 
    
    !MRS.representable && return MRS    
    phi = small_reduce_full_rec(gens(MRS.defining_ideal))
    new_ring = codomain(phi)
    new_ideal = phi(MRS.defining_ideal)
    
    new_ideal_gens = filter!(x->!iszero(x), phi.(gens(MRS.defining_ideal)))
    
    if length(new_ideal_gens) == 0
        new_ideal = ideal(new_ring, [new_ring(0)])
    else
        new_ideal = ideal(new_ideal_gens)
    end

    new_mat = phi.(collect(MRS.representation_matrix))
    new_mat = matrix(new_ring, new_mat)
    new_ineq = unique!(filter(x->!is_unit(x), phi.(MRS.inequations)))
    new_ineq = gens_2_factors(new_ineq)
    
    return MatroidRealizationSpace(new_ideal, new_ineq, new_ring, new_mat, MRS.representable, MRS.F, MRS.char, MRS.q)
    
end



















#####################
# full reduction    #
#####################

# computes the coefficient of v in monomial m. 
function coefficient_monomial(v::RingElem, m::RingElem)
    isone(degree(m,v)) || return "The variable is not a degree 1 factor."
    mf = factor(m)
    mfdict = Dict(mf)
    u = unit(mf)
    not_v = [k^(mfdict[k]) for k in keys(mfdict) if k != v ]
    length(not_v) == 0 ? u : u*prod(not_v)
end

# computes the coefficient of v in f. 
function coefficient_v(v::RingElem, f::RingElem)
    isone(degree(f,v)) || return "degree of variable must be 1"
    withv = [term(f,i) for i in 1:length(f) if v in vars(monomial(f,i))]
    return sum([coefficient_monomial(v,m) for m in withv])
end


function find_solution_v(v::RingElem, Igens::Vector{<:RingElem}, 
                         Sgens::Vector{<:RingElem}, R::MPolyRing) 

    
    with_v_deg_1 = [g for g in Igens if isone(degree(g,v))] 
    length(with_v_deg_1) != 0 || return "can't isolate"

    for f in with_v_deg_1

        den = coefficient_v(v, f)
        fac_den = poly_2_factors(den)
        !issubset(fac_den, Sgens) && continue

        no_v = [term(f,i) for i in 1:length(f) if !(v in vars(monomial(f,i)))]
        iszero(length(no_v)) && continue
              
        h = R(-1)*sum(no_v)
        return h//den
        
    end
    return "can't solve for v"
end


# v is replaced by t in f
function sub_map(v::RingElem, t::RingElem, R::MPolyRing, xs::Vector{<:RingElem}) 
    xs_v = map(x -> x==v ? t : x, xs )    
    return hom(R,FractionField(R), a->a, xs_v)
end



# replace v by t in f, only return the numerator.
function sub_v(v::RingElem, t::RingElem, f::RingElem, R::Ring, xs::Vector{<:RingElem}) 
    m = sub_map(v,t,R,xs) 
    new_f = numerator(m(f))
    return new_f 
end


# removes factors that are in the semigroup generated by Sgens
function clean(f::RingElem, R::MPolyRing, Sgens::Vector{<:RingElem})   
    
    fFactors = factor(f)
    FactorsDict = Dict(fFactors) 
    cleanf_arr = [k^(FactorsDict[k]) for k in keys(FactorsDict) if !(k in Sgens) || is_unit(k)]
    
    length(cleanf_arr) > 0 ? prod(cleanf_arr) : unit(fFactors)
    
end

#variables in ideal
function ideal_vars(Igens::Vector{<:RingElem}) 
    return unique!(vcat([vars(gen) for gen in Igens]...))
end

function n_new_Sgens(x::RingElem, t::RingElem, Sgens::Vector{<:RingElem}, 
                     R::Ring, xs::Vector{<:RingElem}) 
    preSgens = unique!([sub_v(x, t, f, R, xs) for f in Sgens])
    return gens_2_factors(preSgens)
end

function n_new_Igens(x::RingElem, t::RingElem, Igens::Vector{<:RingElem}, 
                     Sgens::Vector{<:RingElem}, R::Ring, xs::Vector{<:RingElem}) 

    preIgens = unique!([clean(sub_v(x, t, f, R, xs), R, Sgens) for f in Igens])
    return filter(x-> x!= R(0), preIgens)
end



function matrix_clear_den_in_col(X::Oscar.MatElem, c::Int)
    Xc = [denominator(f) for f in X[:, c]]
    t = lcm(Xc)
    result = multiply_column!(X, t, c)
    return result

end



function matrix_clear_den(X::Oscar.MatElem)
    rs, cs = size(X)
    for c in 1:cs
        X = matrix_clear_den_in_col(X, c)
    end
    return X
end


function reduce_ideal_one_step(MRS::MatroidRealizationSpace, 
                               elim::Vector{<:RingElem}, 
                               fullyReduced::Bool)

    Igens = gens(MRS.defining_ideal)
    Sgens = MRS.inequations
    R = MRS.ambient_ring
    FR = FractionField(R)
    xs = gens(R)
    X = MRS.representation_matrix
    nr, nc = size(X)
    
    Ivars = ideal_vars(Igens);

    for x in Ivars 
        t = find_solution_v(x, Igens, Sgens, R)
        t isa String && continue
        
        phi = sub_map(x, t, R, xs)
        
	Sgens_new = n_new_Sgens(x, t, Sgens, R, xs);
        Igens_new = n_new_Igens(x, t, Igens, Sgens_new, R, xs);
        push!(elim, x)
        
        phiX = matrix(FR, [phi(X[i,j]) for i in 1:nr, j in 1:nc  ] )
        nX_FR = matrix_clear_den(phiX)
        nX = matrix(R, [numerator(nX_FR[i,j])  for i in 1:nr, j in 1:nc ])
        
        GBnew = collect(groebner_basis(ideal(R, Igens_new)))         
        
        MRS_new = MatroidRealizationSpace(ideal(R, GBnew), Sgens_new, R, nX, MRS.representable, MRS.F, MRS.char, MRS.q )

        
        return (MRS_new, elim, fullyReduced)
    end

    return (MRS, elim, true)

end


function reduce_realization_space(MRS::MatroidRealizationSpace,
                               elim::Vector{RingElem} = Vector{RingElem}(), 
                               fullyReduced::Bool = false) 
    
    
    #If there are no variables left, we don't reduce anything
    if !(typeof(MRS.ambient_ring) <: MPolyRing)
        return MRS
    end

    output = reduce_ideal_one_step(MRS, elim, fullyReduced)
    output isa String && return "Not Realizable 0 in Semigroup"
    (MRS, elim, fullyReduced) = output    
    
    !fullyReduced && return reduce_realization_space(MRS, elim, fullyReduced)
    
    
    R = MRS.ambient_ring
    xs = gens(R)
    cR = coefficient_ring(R)
    X = MRS.representation_matrix
    nr, nc = size(X)
    Igens = gens(MRS.defining_ideal)
    Sgens = MRS.inequations

    
    zero_elim = []        
    for i in 1:length(xs)
        if xs[i] in elim
            push!(zero_elim, 0)
        else
            push!(zero_elim, "x"*string(i) ) 
        end
    end
        
    xnew_str = Vector{String}(filter(x -> x!=0,  zero_elim))    
    
    
    
    if length(xnew_str) == 0
        phi = hom(R, cR, a->a, [cR(0) for i in 1:length(xs)])
        ambR = codomain(phi);
        if length(Igens) == 0
            Inew = ideal(ambR, ambR(0))
        else
            Inew = ideal(ambR, phi.(Igens)); 
        end
        
        if length(Sgens) == 0
            normal_Sgens = Vector{RingElem}[]
        else             
            normal_Sgens = phi.(Sgens)
        end        
        
        
    else
        Rnew, xnew = polynomial_ring(coefficient_ring(R), length(xnew_str)) 
    
        zero_elim_var = []
        j=1
        for i in 1:length(zero_elim)
            if xs[i] in elim
                push!(zero_elim_var, Rnew(0))
            else
                push!(zero_elim_var, xnew[j] ) 
                j+=1
            end
        end
        
        phi = hom(R, Rnew, a->a, zero_elim_var)
        ambR = codomain(phi);
        if length(Igens) == 0
            Inew = ideal(ambR, ambR(0))
        else
            Inew = ideal(ambR, phi.(Igens)); 
        end
        
        if length(Sgens) == 0
            normal_Sgens = Vector{RingElem}[]
        else             
            Sgens_new = phi.(Sgens)
            normal_Sgens = gens_2_factors([normal_form(g, Inew) for g in Sgens_new])
            unique!(normal_Sgens)
        end
    
    end
    
    Xnew = matrix(ambR, [phi(X[i,j]) for i in 1:nr, j in 1:nc])

    #Try to reduce the matrix one last time using the ideal and the inequations
    m, n = size(Xnew)
    A = base_ring(R)
    if Inew!=ideal(A,A(0)) && length(gens(Inew)) > 0
        for i in 1:m
            for j in 1:n
                if A == ZZ
                    Xnew[i,j] = mod(Xnew[i,j], gens(Inew)[1])
                else
                    Xnew[i,j] = reduce(Xnew[i,j], gens(Inew))
                end
            end
        end
    end

    for j in 1:n
        g = gcd(Xnew[:,j]...)
        factors = poly_2_factors(g)
        for f in factors
            if f in normal_Sgens
                for i in 1:m
                    Xnew[i,j] = Xnew[i,j]/f
                end
            end
        end
    end

    MRS_new = MatroidRealizationSpace(Inew, normal_Sgens, ambR, Xnew, MRS.representable, MRS.F, MRS.char, MRS.q)

    return MRS_new
end




























































############
# Old      #
############



#function reduce_ideal_one_step(Igens::Vector{<:RingElem}, Sgens::Vector{<:RingElem}, R::Ring, 
#                               xs::Vector{<:RingElem}, elim::Vector{<:RingElem}, fullyReduced::Bool)
    
#    if R(0) in Sgens
#        return "Not realizable 0 in Semigroup"
#    else
#        Ivars = ideal_vars(Igens); 
        
        
#        for x in Ivars 
#           t = find_solution_v(x, Igens, Sgens, R)
#            if t isa String
#                continue
#            else 
#                Sgens_new = n_new_Sgens(x, t, Sgens, R, xs);
#                Igens_new = n_new_Igens(x, t, Igens, Sgens_new, R, xs);
#                push!(elim, x)
#                return (Igens_new, Sgens_new, R, xs, elim, fullyReduced)
#            end
   
#        end
#    return (Igens, Sgens, R, xs, elim, true)
#    end 
#end


#function reduce_ideal_full_rec(Igens::Vector{<:RingElem}, Sgens::Vector{<:RingElem}, R::MPolyRing, 
#                           xs::Vector{<:RingElem}, elim::Vector{RingElem} = Vector{RingElem}(), 
#                           fullyReduced::Bool = false) 
    
    
#    output = reduce_ideal_one_step(Igens, Sgens, R, xs, elim, fullyReduced)
#    output isa String && return "Not Realizable 0 in Semigroup"
#    (Igens, Sgens, R, xs, elim, fullyReduced) = output    
    
#    !fullyReduced && return reduce_ideal_full_rec(Igens, Sgens, R, xs, elim, fullyReduced)
    
#    zero_elim = []        
#    for i in 1:length(xs)
#        if xs[i] in elim
#            push!(zero_elim, 0)
#        else
#            push!(zero_elim, "x"*string(i) ) 
#        end
#    end
        
#    xnew_str = Vector{String}(filter(x -> x!=0,  zero_elim))    
    
#    cR = coefficient_ring(R)
    
#    if length(xnew_str) == 0
        
#        phi = hom(R, cR, a->a, [cR(0) for i in 1:length(xs)])
#        return (phi.(Igens), phi.(Sgens), cR)
   
#    end
    

#    Rnew, xnew = PolynomialRing(coefficient_ring(R), xnew_str) 
    
#    zero_elim_var = []
#    j=1
#    for i in 1:length(zero_elim)
#        if xs[i] in elim
#            push!(zero_elim_var, Rnew(0))
#        else
#            push!(zero_elim_var, xnew[j] ) 
#            j+=1
#        end
#    end
        
#    phi = hom(R, Rnew, a->a, zero_elim_var)
    
#    return (phi.(Igens), phi.(Sgens), Rnew)
#end

#function reduce_ideal_full(MRS::MatroidRealizationSpace)
#    R= MRS.ambient_ring
#    xs=gens(R)
#    return reduce_ideal_full_rec(gens(MRS.defining_ideal), MRS.inequations, R, xs)

#end


