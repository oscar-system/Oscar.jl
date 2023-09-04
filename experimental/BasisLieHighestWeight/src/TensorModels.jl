using Oscar

function kron(A::SMat{ZZRingElem}, B::SMat{ZZRingElem})::SMat{ZZRingElem}
    """
    Computes the Kronecker-product of A and B
    """
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
    return res
end

# temprary fix sparse in Oscar does not work
function tensorProduct(A::SMat{ZZRingElem}, B::SMat{ZZRingElem})::SMat{ZZRingElem}
    temp_mat = kron(A, spid(sz(B))) + kron(spid(sz(A)), B)
    res = sparse_matrix(ZZ, nrows(A)*nrows(B), ncols(A)*ncols(B))
    for i in 1:nrows(temp_mat)
        setindex!(res, getindex(temp_mat, i), i)
    end
    return res
end


spid(n::Int) = identity_matrix(SMat, ZZ, n)::SMat{ZZRingElem}
sz(A::SMat{ZZRingElem}) = nrows(A)::Int #size(A)[1]
#tensorProduct(A, B) = kron(A, spid(sz(B))) + kron(spid(sz(A)), B)
tensorProducts(As, Bs) = (AB->tensorProduct(AB[1], AB[2])).(zip(As, Bs))
tensorPower(A, n) = (n == 1) ? A : tensorProduct(tensorPower(A, n-1), A)
tensorPowers(As, n) = (A->tensorPower(A, n)).(As)

function tensorMatricesForOperators(lie_algebra::GAP.Obj, highest_weight::Vector{ZZRingElem}, 
                                    operators::GAP.Obj)::Vector{SMat{ZZRingElem}}
    """
    Calculates the matrices g_i corresponding to the operator ops[i].
    """
    matrices_of_operators = []
    for i in 1:length(highest_weight)
        if highest_weight[i] <= 0
            continue
        end
        wi = convert(Vector{ZZRingElem}, Int.(1:length(highest_weight) .== i)) # i-th fundamental weight
        _matrices_of_operators = matricesForOperators(lie_algebra, wi, operators)
        _matrices_of_operators = tensorPowers(_matrices_of_operators, highest_weight[i])
        matrices_of_operators = matrices_of_operators == [] ? _matrices_of_operators : 
                                  tensorProducts(matrices_of_operators, _matrices_of_operators)
    end
    return matrices_of_operators
end
