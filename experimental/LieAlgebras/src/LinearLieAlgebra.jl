@attributes mutable struct LinearLieAlgebra{C<:RingElement} <: LieAlgebra{C}
  R::Ring
  n::Int  # the n of the gl_n this embeds into
  dim::Int
  basis::Vector{MatElem{C}}
  s::Vector{Symbol}

  function LinearLieAlgebra{C}(
    R::Ring, n::Int, basis::Vector{<:MatElem{C}}, s::Vector{Symbol}; cached::Bool=true
  ) where {C<:RingElement}
    return get_cached!(
      LinearLieAlgebraDict, (R, n, basis, s), cached
    ) do
      @req all(b -> size(b) == (n, n), basis) "Invalid basis element dimensions."
      @req length(s) == length(basis) "Invalid number of basis element names."
      new{C}(R, n, length(basis), basis, s)
    end::LinearLieAlgebra{C}
  end
end

const LinearLieAlgebraDict = CacheDictType{
  Tuple{Ring,Int,Vector{<:MatElem},Vector{Symbol}},LinearLieAlgebra
}()

struct LinearLieAlgebraElem{C<:RingElement} <: LieAlgebraElem{C}
  parent::LinearLieAlgebra{C}
  mat::MatElem{C}
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LinearLieAlgebraElem{C}}) where {C<:RingElement} = LinearLieAlgebra{C}

elem_type(::Type{LinearLieAlgebra{C}}) where {C<:RingElement} = LinearLieAlgebraElem{C}

parent(x::LinearLieAlgebraElem{C}) where {C<:RingElement} = x.parent

base_ring(L::LinearLieAlgebra{C}) where {C<:RingElement} = L.R::parent_type(C)

dim(L::LinearLieAlgebra{C}) where {C<:RingElement} = L.dim

function matrix_repr_basis(L::LinearLieAlgebra{C}) where {C<:RingElement}
  return Vector{dense_matrix_type(C)}(L.basis)
end

function matrix_repr_basis(L::LinearLieAlgebra{C}, i::Int) where {C<:RingElement}
  return (L.basis[i])::dense_matrix_type(C)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::LinearLieAlgebra{C}) where {C<:RingElement}
  print(io, "LinearLieAlgebra (⊆ gl_$(V.n)) over ")
  print(IOContext(io, :compact => true), base_ring(V))
end

function symbols(L::LinearLieAlgebra{C}) where {C<:RingElement}
  return L.s
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (L::LinearLieAlgebra{C})(m::MatElem{C}) where {C<:RingElement}
  if size(m) == (L.n, L.n)
    m = coefficient_vector(m, matrix_repr_basis(L))
  end
  @req size(m) == (1, dim(L)) "Invalid matrix dimensions."
  return elem_type(L)(L, m)
end

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function Generic.matrix_repr(x::LinearLieAlgebraElem{C}) where {C<:RingElement}
  return sum(c * b for (c, b) in zip(x.mat, matrix_repr_basis(parent(x))))
end

function bracket(
  x::LinearLieAlgebraElem{C}, y::LinearLieAlgebraElem{C}
) where {C<:RingElement}
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

function lie_algebra(
  R::Ring, n::Int, basis::Vector{<:MatElem{C}}, s::Vector{<:VarName}; cached::Bool=true
) where {C<:RingElement}
  return LinearLieAlgebra{elem_type(R)}(R, n, basis, Symbol.(s); cached)
end

function general_linear_lie_algebra(R::Ring, n::Int)
  basis = [(b = zero_matrix(R, n, n); b[i, j] = 1; b) for i in 1:n for j in 1:n]
  s = ["x_$(i)_$(j)" for i in 1:n for j in 1:n]
  L = lie_algebra(R, n, basis, s)
  set_attribute!(L, :type => :general_linear)
  return L
end

function special_linear_lie_algebra(R::Ring, n::Int)
  basis_e = [(b = zero_matrix(R, n, n); b[i, j] = 1; b) for i in 1:n for j in (i + 1):n]
  basis_f = [(b = zero_matrix(R, n, n); b[j, i] = 1; b) for i in 1:n for j in (i + 1):n]
  basis_h = [
    (b = zero_matrix(R, n, n); b[i, i] = 1; b[i + 1, i + 1] = -1; b) for i in 1:(n - 1)
  ]
  s_e = ["e_$(i)_$(j)" for i in 1:n for j in (i + 1):n]
  s_f = ["f_$(i)_$(j)" for i in 1:n for j in (i + 1):n]
  s_h = ["h_$(i)" for i in 1:(n - 1)]
  L = lie_algebra(R, n, [basis_e; basis_f; basis_h], [s_e; s_f; s_h])
  set_attribute!(L, :type => :special_linear)
  return L
end

function special_orthogonal_lie_algebra(R::Ring, n::Int)
  basis = [
    (b = zero_matrix(R, n, n); b[i, j] = 1; b[j, i] = -1; b) for i in 1:n for j in (i + 1):n
  ]
  s = ["x_$(i)_$(j)" for i in 1:n for j in (i + 1):n]
  L = lie_algebra(R, n, basis, s)
  set_attribute!(L, :type => :special_orthogonal)
  return L
end
