#################################################
#
# Abstract parent type
#
#################################################

abstract type LieAlgebra{C <: RingElement} <: FPModule{C} end

abstract type LieAlgebraElem{C <: RingElement} <: FPModuleElem{C} end

# To be implemented by subtypes:
# parent_type(::Type{MyLieAlgebraElem{C}})
# elem_type(::Type{MyLieAlgebra{C}})
# parent(x::MyLieAlgebraElem{C})
# base_ring(L::MyLieAlgebra{C})
# dim(L::MyLieAlgebra{C})
# Generic._matrix(x::MyLieAlgebraElem{C})
# Base.show(io::IO, x::MyLieAlgebra{C})
# Base.show(io::IO, x::MyLieAlgebraElem{C})
# symbols(L::MyLieAlgebra{C})
# bracket(x::MyLieAlgebraElem{C}, y::MyLieAlgebraElem{C})


###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_ring(x::LieAlgebraElem{C}) where {C <: RingElement} = base_ring(parent(x))

ngens(L::LieAlgebra{C}) where {C <: RingElement} = dim(L)

gens(L::LieAlgebra{C}) where {C <: RingElement} = [gen(L, i)::elem_type(L) for i in 1:dim(L)]

function gen(L::LieAlgebra{C}, i::Int) where {C <: RingElement}
    R = base_ring(L)
    return L([(j == i ? one(R) : zero(R)) for j in 1:dim(L)])
end

function Generic.rels(_::LieAlgebra{C}) where {C <: RingElement}
    # there are no relations in a vector space
    return Vector{dense_matrix_type(C)}(undef, 0)
end


###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (L::LieAlgebra{C})() where {C <: RingElement}
    mat = zero_matrix(base_ring(L), 1, dim(L))
    return elem_type(L)(L, mat)
end

function (L::LieAlgebra{C})(v::Vector{Int}) where {C <: RingElement}
    return L(base_ring(L).(v))
end

function (L::LieAlgebra{C})(v::Vector{C}) where {C <: RingElement}
    length(v) == dim(L) || error("Length of vector does not match number of generators.")
    mat = matrix(base_ring(L), 1, length(v), v)
    return elem_type(L)(L, mat)
end

function (L::LieAlgebra{C})(mat::MatElem{C}) where {C <: RingElement}
    size(mat) == (1, dim(L)) || error("Invalid matrix dimensions.")
    return elem_type(L)(L, mat)
end

function (L::LieAlgebra{C})(v::SRow{C}) where {C <: RingElement}
    mat = dense_row(v, dim(L))
    return elem_type(L)(L, mat)
end

function (L::LieAlgebra{C})(v::LieAlgebraElem{C}) where {C <: RingElement}
    L == parent(v) || error("Incompatible modules.")
    return v
end


###############################################################################
#
#   Arithmetic operations
#
###############################################################################

# Vector space operations get inherited from FPModule


###############################################################################
#
#   Comparison functions
#
###############################################################################

# Overwrite the equality of FPModule to be used for CacheDicts
function Base.:(==)(L1::LieAlgebra{C}, L2::LieAlgebra{C}) where {C <: RingElement}
    return L1 === L2
end


###############################################################################
#
#   Constructor
#
###############################################################################

function general_linear_liealgebra(R::Ring, n::Int)
    basis = [(b = zero_matrix(R, n, n); b[i, j] = 1; b) for i in 1:n for j in 1:n]
    s = ["x_$(i)_$(j)" for i in 1:n for j in 1:n]
    L = liealgebra(R, n, basis, s)
    set_attribute!(L, :type, :general_linear)
    return L
end

function special_linear_liealgebra(R::Ring, n::Int)
    basis_e = [(b = zero_matrix(R, n, n); b[i, j] = 1; b) for i in 1:n for j in i+1:n]
    basis_f = [(b = zero_matrix(R, n, n); b[j, i] = 1; b) for i in 1:n for j in i+1:n]
    basis_h = [(b = zero_matrix(R, n, n); b[i, i] = 1; b[i+1, i+1] = -1; b) for i in 1:n-1]
    s_e = ["e_$(i)_$(j)" for i in 1:n for j in i+1:n]
    s_f = ["f_$(i)_$(j)" for i in 1:n for j in i+1:n]
    s_h = ["h_$(i)" for i in 1:n-1]
    L = liealgebra(R, n, [basis_e; basis_f; basis_h], [s_e; s_f; s_h])
    set_attribute!(L, :type, :special_linear)
    return L
end

function special_orthogonal_liealgebra(R::Ring, n::Int)
    basis = [(b = zero_matrix(R, n, n); b[i, j] = 1; b[j, i] = -1; b) for i in 1:n for j in i+1:n]
    s = ["x_$(i)_$(j)" for i in 1:n for j in i+1:n]
    L = liealgebra(R, n, basis, s)
    set_attribute!(L, :type, :special_orthogonal)
    return L
end
