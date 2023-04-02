@attributes mutable struct LinearLieAlgebra{C <: RingElement} <: LieAlgebra{C}
    R::Ring
    n::Int  # the n of the gl_n this embeds into
    dim::Int
    basis::Vector{MatElem{C}}
    s::Vector{Symbol}

    function LinearLieAlgebra{C}(
        R::Ring,
        n::Int,
        basis::Vector{<:MatElem{C}},
        s::Vector{Symbol};
        cached::Bool=true,
    ) where {C <: RingElement}
        return get_cached!(LinearLieAlgebraDict, (R, n, basis, s), cached) do
            all(b -> size(b) == (n, n), basis) || error("Invalid basis element dimensions.")
            length(s) == length(basis) || error("Invalid number of basis element names.")
            new{C}(R, n, length(basis), basis, s)
        end::LinearLieAlgebra{C}
    end
end

const LinearLieAlgebraDict = CacheDictType{Tuple{Ring, Int, Vector{<:MatElem}, Vector{Symbol}}, LinearLieAlgebra}()

struct LinearLieAlgebraElem{C <: RingElement} <: LieAlgebraElem{C}
    parent::LinearLieAlgebra{C}
    mat::MatElem{C}
end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LinearLieAlgebraElem{C}}) where {C <: RingElement} = LinearLieAlgebra{C}

elem_type(::Type{LinearLieAlgebra{C}}) where {C <: RingElement} = LinearLieAlgebraElem{C}

parent(x::LinearLieAlgebraElem{C}) where {C <: RingElement} = x.parent

base_ring(L::LinearLieAlgebra{C}) where {C <: RingElement} = L.R

dim(L::LinearLieAlgebra{C}) where {C <: RingElement} = L.dim

@inline function Generic._matrix(x::LinearLieAlgebraElem{C}) where {C <: RingElement}
    return (x.mat)::dense_matrix_type(C)
end

function matrix_repr_basis(L::LinearLieAlgebra{C}) where {C <: RingElement}
    return Vector{dense_matrix_type(C)}(L.basis)
end

function matrix_repr_basis(L::LinearLieAlgebra{C}, i::Int) where {C <: RingElement}
    return (L.basis[i])::dense_matrix_type(C)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::LinearLieAlgebra{C}) where {C <: RingElement}
    print(io, "LinearLieAlgebra (⊆ gl_$(V.n)) over ")
    print(IOContext(io, :compact => true), base_ring(V))
end

function Base.show(io::IO, x::LinearLieAlgebraElem{C}) where {C <: RingElement}
    Base.show(io, matrix_repr(x))
end

function symbols(L::LinearLieAlgebra{C}) where {C <: RingElement}
    return L.s
end


###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (L::LinearLieAlgebra{C})(m::MatElem{C}) where {C <: RingElement}
    if size(m) == (L.n, L.n)
        m = coefficient_vector(m, matrix_repr_basis(L))
    end
    size(m) == (1, dim(L)) || error("Invalid matrix dimensions.")
    return elem_type(L)(L, m)
end


###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function Generic.matrix_repr(x::LinearLieAlgebraElem{C}) where {C <: RingElement}
    return sum(c * b for (c, b) in zip(x.mat, matrix_repr_basis(parent(x))))
end

function bracket(x::LinearLieAlgebraElem{C}, y::LinearLieAlgebraElem{C}) where {C <: RingElement}
    check_parent(x, y)
    L = parent(x)
    x_mat = matrix_repr(x)
    y_mat = matrix_repr(y)
    return L(x_mat * y_mat - y_mat * x_mat)
end


###############################################################################
#
#   Constructor
#
###############################################################################

function liealgebra(
    R::Ring,
    n::Int,
    basis::Vector{<:MatElem{C}},
    s::Vector{<:Union{AbstractString, Char, Symbol}};
    cached::Bool=true,
) where {C <: RingElement}
    return LinearLieAlgebra{elem_type(R)}(R, n, basis, Symbol.(s); cached)
end
