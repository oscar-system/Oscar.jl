export ComplexReflection
export scalar_product
export unitary_reflection
export hyperplane
export coroot
export eigenvalue

struct ComplexReflection{T <: QQAlgFieldElem}
    base_ring::QQAlgField #this should be the parent of T
    matrix::MatElem{T}
    root::AbstractAlgebra.Generic.FreeModuleElem{T}
    coroot::AbstractAlgebra.Generic.FreeModuleElem{T}
    hyperplane::MatElem{T}
    eigenvalue::T
    order::Int
    root_norm_squared::T
    is_unitary::Bool
end

# Printing
function Base.show(io::IO, ::MIME"text/plain", w::ComplexReflection)
    print(io, "Complex reflection of order ", w.order, " with root ", w.root)
end

# Getter functions
function base_ring(w::ComplexReflection)
    return w.base_ring
end

function matrix(w::ComplexReflection)
    return w.matrix
end

function root(w::ComplexReflection)
    return w.root
end

function coroot(w::ComplexReflection)
    return w.coroot
end

function hyperplane(w::ComplexReflection)
    return w.hyperplane
end

function eigenvalue(w::ComplexReflection)
    return w.eigenvalue
end

function order(w::ComplexReflection)
    return w.eigenvalue
end

function is_unitary(w::ComplexReflection)
    return w.is_unitrary
end

function scalar_product(v::AbstractAlgebra.Generic.FreeModuleElem{T}, w::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: QQAlgFieldElem

    V = parent(v)
    K = base_ring(V)
    n = dim(V)
    s = zero(K)
    for i=1:n
        s += v[i]*w[i]
    end
    return s
end

function unitary_reflection(root::AbstractAlgebra.Generic.FreeModuleElem{T}, zeta::T, order::Int) where T <: QQAlgFieldElem

    V = parent(root)
    K = base_ring(V)
    n = dim(V)

    root_norm_squared = scalar_product(root,root)
    c = (1-zeta)//root_norm_squared
    coroot = c*root

    basis = gens(V)

    w = matrix([ basis[i] - coroot * root[i] for i=1:n ])

    I = identity_matrix(K, n)

    hyp = kernel(I-w)
    
    return ComplexReflection{typeof(K(1))}(K, w, root, coroot, hyp, zeta, order, root_norm_squared, true)

end

function unitary_reflection(root::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: QQAlgFieldElem

    V = parent(root)
    K = base_ring(V)

    zeta = K(-1)
    
    return unitary_reflection(root, zeta, 2)

end