
export is_unitary

function is_unitary(M::AbstractAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem})

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

function is_unitary(M::AbstractAlgebra.Generic.MatSpaceElem{AbsNonSimpleNumFieldElem})

    if !is_square(M)
        return false
    end

    # convert base field to absolute and convert matrix as well
    K = base_ring(M)
    n = ncols(M)
    L, f = absolute_simple_field(K)
    ML = matrix(L, n, n, [f(x) for x in M])

    return is_unitary(ML)

end