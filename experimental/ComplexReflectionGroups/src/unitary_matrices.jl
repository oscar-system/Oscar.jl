export is_unitary
export is_orthogonal

# Unitary matrices
#
# Ulrich Thiel, 2023

function is_orthogonal(M::MatElem)
    if !is_square(M)
        return false
    end

    K = base_ring(M)
    n = ncols(M)

    return M*transpose(M) == identity_matrix(K,n)
end

function is_unitary(M::QQMatrix)
    
    return is_orthogonal(M)

end

function complex_conjugation(K::QQField)

    return identity_map(K)

end

function complex_conjugation(K::NumField)

    L, f = absolute_simple_field(K)
    g = inv(f)
    conj = complex_conjugation(L)
    return compose(compose(g, conj), f)

end

function is_unitary(M::MatElem{T}) where T <: NumFieldElem

    if !is_square(M)
        return false
    end

    # create the conjugate transpose of M
    K = base_ring(M)
    n = ncols(M)
    conj = complex_conjugation(K)
    Mct = transpose(matrix(K, n, n, [conj(x) for x in M]))

    return M*Mct == identity_matrix(K,n)

end

function is_unitary(M::MatrixGroupElem{T}) where T <: QQAlgFieldElem
    return is_unitary(matrix(M))
end

function is_unitary(G::MatrixGroup{T}) where T <: QQAlgFieldElem
    for g in gens(G)
        if is_unitary(g) == false
            return false
        end
    end
    return true
end