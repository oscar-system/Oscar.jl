export is_unitary

# Unitary matrices
#
# Ulrich Thiel, 2023

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
    g = inv(f)
    ML = matrix(L, n, n, [g(x) for x in M])

    return is_unitary(ML)

end

function is_unitary(M::MatrixGroupElem{AbsSimpleNumFieldElem})
    return is_unitary(matrix(M))
end

function is_unitary(G::MatrixGroup{AbsSimpleNumFieldElem})
    for g in gens(G)
        if is_unitary(g) == false
            return false
        end
    end
    return true
end