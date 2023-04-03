# ?

using Oscar
#using SparseArrays

# TODO: make the first one a symmetric product, or reduce more generally

function kron(A, B)
    res = sparse_matrix(ZZ, nrows(A)*nrows(B), ncols(A)*ncols(B))
    for i in 1:nrows(B)
        for j in 1:nrows(A)
            new_row_tuples = Vector{Tuple{Int, ZZRingElem}}([(1,ZZ(0))])
            for (index_A, element_A) in union(getindex(A, j))
                for (index_B, element_B) in union(getindex(B, i))
                    push!(new_row_tuples, ((index_A-1)*ncols(B)+index_B, element_A*element_B))
                end
            end
            new_row = sparse_row(ZZ, new_row_tuples)
            setindex!(res, new_row, (j-1)*nrows(B)+i)
        end
    end
    #println("ncols(res): ", ncols(res))
    #println("nrows(res): ", nrows(res))
    
    return res
end

# temprary fix sparse in Oscar does not work
function tensorProduct(A, B)
    temp_mat = kron(A, spid(sz(B))) + kron(spid(sz(A)), B)
    res = sparse_matrix(ZZ, nrows(A)*nrows(B), ncols(A)*ncols(B))
    for i in 1:nrows(temp_mat)
        setindex!(res, getindex(temp_mat, i), i)
    end
    return res
end


spid(n) = identity_matrix(SMat, ZZ, n)
sz(A) = nrows(A) #size(A)[1]
#tensorProduct(A, B) = kron(A, spid(sz(B))) + kron(spid(sz(A)), B)
tensorProducts(As, Bs) = (AB->tensorProduct(AB[1], AB[2])).(zip(As, Bs))
tensorPower(A, n) = (n == 1) ? A : tensorProduct(tensorPower(A, n-1), A)
tensorPowers(As, n) = (A->tensorPower(A, n)).(As)


function tensorMatricesForOperators(L, hw, ops)
    """
    Calculates the matrices g_i corresponding to the operator ops[i].
    """
    #println("hw: ", hw)
    mats = []

    for i in 1:length(hw)
        #println("hw[i]: ", hw[i])
        if hw[i] <= 0
            continue
        end
        wi = Int.(1:length(hw) .== i) # i-th fundamental weight
        _mats = matricesForOperators(L, wi, ops)
        _mats = tensorPowers(_mats, hw[i])
        mats = mats == [] ? _mats : tensorProducts(mats, _mats)
        #println(spdiagm(0 => [ZZ(1) for _ in 1:5])agm)
        #display(mats)
    end
    return mats
end