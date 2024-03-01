export ComplexReflection
export scalar_product
export unitary_reflection
export hyperplane
export eigenvalue
export is_complex_reflection_with_data
export is_complex_reflection
export is_root_of_unity_with_data
export is_root_of_unity
export is_complex_reflection_group
export complex_reflections

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
    return w.order
end

function is_unitary(w::ComplexReflection)
    return w.is_unitary
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

function is_root_of_unity_with_data(x::QQFieldElem)

    if x == 0
        return false, 0
    elseif x == 1
        return true, 1
    elseif x == -1
        return true, 2
    else
        return false, 0
    end

end

function is_root_of_unity_with_data(x::QQAlgFieldElem)

    if x == 0
        return false, 0
    elseif x == 1
        return true, 1
    elseif x == -1
        return true, 2
    else
        p = absolute_minpoly(x)
        if !is_cyclotomic_polynomial(p)
            return false, 0
        end
        candidates = euler_phi_inv(degree(p))
        for n in sort(candidates)
            if x^n == 1
                return true, n
            end
        end
    end
end

function is_root_of_unity(x::QQAlgFieldElem)
    b,n = is_root_of_unity_with_data(x)
    return b
end

function is_complex_reflection_with_data(w::MatElem{T}) where T <: QQAlgFieldElem
    
    if !is_square(w)
        return false, nothing
    end

    n = number_of_rows(w)
    K = base_ring(w)
    I = identity_matrix(K, n)
    V = vector_space(K,n)
    f = module_homomorphism(V,V,I-w)

    H, Hincl = kernel(f)

    if dim(H) != n-1
        return false, nothing
    end

    R, Rincl = image(f)

    # If H \cap R != {0} then w is not of finite order: 
    # if w is of finite order, then image(1-w) is the orthogonal complement of ker(1-w)
    # with respect to a w-invariant scalar product, which exists since w is of finite order.
    # I am not sure if this is an if and only if, but we will determine the order of w 
    # below anyways.
    HcapR, HcapRincl = intersect(H,R)

    if dim(HcapR) > 0
        return false, nothing
    end

    # Now, if H \cap R = {0}, then R must be 1-dimensional and we take as root a (the) basis
    # vector
    alpha = Rincl(gens(R)[1])

    # Now, we need to determine the scalar zeta by which w acts on alpha (we do not yet know
    # whether w is a root of unity (which is equivalent to w being of finite order).
    alpha_w = alpha*w

    # We need a non-zero entry of alpha for this
    i=1
    while alpha[i] == 0
        i += 1
    end
    zeta = alpha_w[i]//alpha[i]

    # if zeta = 0, then w is not invertible
    # this happens for example for w = matrix(QQ,2,2,[1 0; 0 0])
    if zeta == 0
        return false, nothing
    end

    # zeta needs to be of finite order (so that w is of finite order), i.e. zeta needs
    # to be a root of unity
    b, d = is_root_of_unity_with_data(zeta)

    if b == false
        return false, nothing
    end

    # We take the normalized coroot 
    alpha_norm_squared = scalar_product(alpha,alpha)
    c = (1-zeta)//alpha_norm_squared
    alpha_check = c*alpha

    w_data = ComplexReflection(K, w, alpha, alpha_check, matrix(Hincl.(gens(H))), zeta, d, alpha_norm_squared, is_unitary(w))

    return true, w_data

end

function is_complex_reflection_with_data(w::MatrixGroupElem{T}) where T <: QQAlgFieldElem
    return is_complex_reflection_with_data(matrix(w))
end

function is_complex_reflection(w::MatElem{T}) where T <: QQAlgFieldElem

    b,data = is_complex_reflection_with_data(w)
    return b

end

function is_complex_reflection(w::MatrixGroupElem{T}) where T <: QQAlgFieldElem
    return is_complex_reflection(matrix(w))
end


function is_complex_reflection_group(G::MatrixGroup{T}) where T <: QQAlgFieldElem
    if has_attribute(G, :is_complex_reflection_group)
        return get_attribute(G, :is_complex_reflection_group)
    end
    
    for g in gens(G)
        if !is_complex_reflection(g)
            set_attribute!(G, :is_complex_reflection_group, false)
            return false
        end
    end

    set_attribute!(G, :is_complex_reflection_group, true)
    return true
end

function complex_reflections(G::MatrixGroup{T}) where T <: QQAlgFieldElem
   
    if has_attribute(G, :complex_reflections)
        return get_attribute(G, :complex_reflections)
    end

    refls = Set(ComplexReflection[])

    # This is not efficient yet: we should loop only over conjugacy classes
    for g in G
        b,g_data = is_complex_reflection_with_data(g)
        if b 
            push!(refls, g_data)
        end
    end

    set_attribute!(G, :complex_reflections, refls)

    return refls

end