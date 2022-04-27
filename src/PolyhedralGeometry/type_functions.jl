################################################################################
######## Scalar types
################################################################################


function detect_scalar_type(n::Type{T}, p::Polymake.BigObject) where T<:Union{Polyhedron, Cone, PolyhedralFan, SubdivisionOfPoints, PolyhedralComplex}
    scalar_regexp = match(r"[^<]*<(.*)>[^>]*", String(Polymake.type_name(p)))
    typename = scalar_regexp[1]
    return scalar_type_to_oscar[typename]
end

scalar_type(::Union{Polyhedron{T}, Cone{T}, Hyperplane{T}, Halfspace{T}}) where T<:scalar_types = T

################################################################################
######## SubObjectIterator
################################################################################

# Matrices with rational only elements
for (sym, name) in (("linear_inequality_matrix", "Linear Inequality Matrix"), ("affine_inequality_matrix", "Affine Inequality Matrix"), ("linear_equation_matrix", "Linear Equation Matrix"), ("affine_equation_matrix", "Affine Equation Matrix"))
    M = Symbol(sym)
    _M = Symbol(string("_", sym))
    @eval begin
        $M(iter::SubObjectIterator{<:Union{Halfspace{fmpq}, Hyperplane{fmpq}, Polyhedron{fmpq}, Cone{fmpq}, Pair{Matrix{fmpq}, fmpq}}}) = matrix(QQ, Matrix{fmpq}($_M(Val(iter.Acc), iter.Obj; iter.options...)))
        $M(iter::SubObjectIterator{<:Union{Halfspace{T}, Hyperplane{T}, Polyhedron{T}, Cone{T}, Pair{Matrix{T}, T}}}) where T<:scalar_types = Matrix{T}($_M(Val(iter.Acc), iter.Obj; iter.options...))
        $M(iter::SubObjectIterator{<:Union{Halfspace{nf_elem}, Hyperplane{nf_elem}, Polyhedron{nf_elem}, Cone{nf_elem}}}) = Matrix{nf_scalar}($_M(Val(iter.Acc), iter.Obj; iter.options...))
        $_M(::Any, ::Polymake.BigObject) = throw(ArgumentError(string($name, " not defined in this context.")))
    end
end

function halfspace_matrix_pair(iter::SubObjectIterator{<:Union{Halfspace{fmpq}, Hyperplane{fmpq}, Polyhedron{fmpq}, Cone{fmpq}, Pair{Matrix{fmpq}, fmpq}}})
    try
        h = affine_matrix_for_polymake(iter)
        return (A = matrix(QQ, Matrix{fmpq}(h[:, 2:end])), b = -h[:, 1])
    catch e
        throw(ArgumentError("Halfspace-Matrix-Pair not defined in this context."))
    end
end

function halfspace_matrix_pair(iter::SubObjectIterator{<:Union{Halfspace{T}, Hyperplane{T}, Polyhedron{T}, Cone{T}, Pair{Matrix{T}, T}}}) where T<:scalar_types
    try
        h = affine_matrix_for_polymake(iter)
        return (A = Matrix{T}(h[:, 2:end]), b = Vector{T}(-h[:, 1]))
    catch e
        throw(ArgumentError("Halfspace-Matrix-Pair not defined in this context."))
    end
end

function halfspace_matrix_pair(iter::SubObjectIterator{<:Union{Halfspace{nf_elem}, Hyperplane{nf_elem}, Polyhedron{nf_elem}, Cone{nf_elem}}})
    try
        h = affine_matrix_for_polymake(iter)
        return (A = Matrix{nf_scalar}(h[:, 2:end]), b = Vector{nf_scalar}(-h[:, 1]))
    catch e
        throw(ArgumentError("Halfspace-Matrix-Pair not defined in this context."))
    end
end
