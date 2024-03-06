# This file implements some basics to work with the Hermitian inner product over
# number fields.

# Ulrich Thiel, 2024

export is_unitary
export is_orthogonal
export complex_conjugation2 


###########################################################################################
# Complex Conjugation 
###########################################################################################
# The key component is complex conjugation. We extend the existing function 
# complex_conjugation(K::AbsSimpleNumField) to more cases of number fields. 
# The function does unfortunately only work for totally complex fields: if the field is 
# totally real, it could simply return the identity. So, we add a new function 
# complex_conjugation2 for the moment.
function complex_conjugation2(K::QQField)

    return identity_map(K)

end

function complex_conjugation2(K::AbsSimpleNumField)
    if is_totally_real(K)
        return identity_map(K)
    else
        return complex_conjugation(K)
    end
end

# There is a function for simple number fields. For a genereal number field we convert
# to an absolute simple and pull conjugation through the isomorphism.
function complex_conjugation2(K::NumField)

    L, f = absolute_simple_field(K)
    g = inv(f)
    conj = complex_conjugation2(L)
    return compose(compose(g, conj), f)

end


###########################################################################################
# Hermitian scalar product
###########################################################################################
function scalar_product(v::AbstractAlgebra.Generic.FreeModuleElem{T}, w::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: QQAlgFieldElem

    V = parent(v)
    K = base_ring(V)
    n = dim(V)
    s = zero(K)
    conj = complex_conjugation2(K)
    for i=1:n
        s += v[i]*conj(w[i])
    end
    return s
end

###########################################################################################
# Complex conjugate of a vector
###########################################################################################
function complex_conjugation(v::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: QQAlgFieldElem

    V = parent(v)
    K = base_ring(V)
    n = dim(V)
    conj = complex_conjugation2(K)
    return V([conj(v[i]) for i=1:n])

end

###########################################################################################
# Unitary matrices
###########################################################################################
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

function is_unitary(M::MatElem{T}) where T <: NumFieldElem

    if !is_square(M)
        return false
    end

    # create the conjugate transpose of M
    K = base_ring(M)
    n = ncols(M)
    conj = complex_conjugation2(K)
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