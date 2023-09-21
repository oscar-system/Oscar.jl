@attributes mutable struct LinearLieAlgebra{C<:RingElement} <: LieAlgebra{C}
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
    check::Bool=true,
  ) where {C<:RingElement}
    return get_cached!(
      LinearLieAlgebraDict, (R, n, basis, s), cached
    ) do
      @req all(b -> size(b) == (n, n), basis) "Invalid basis element dimensions."
      @req length(s) == length(basis) "Invalid number of basis element names."
      L = new{C}(R, n, length(basis), basis, s)
      if check
        @req all(b -> all(e -> parent(e) === R, b), basis) "Invalid matrices."
        # TODO: make work
        # for xi in basis(L), xj in basis(L)
        #   @req (xi * xj) in L
        # end
      end
      return L
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

parent(x::LinearLieAlgebraElem) = x.parent

coefficient_ring(L::LinearLieAlgebra{C}) where {C<:RingElement} = L.R::parent_type(C)

dim(L::LinearLieAlgebra) = L.dim

@doc raw"""
    matrix_repr_basis(L::LinearLieAlgebra{C}) -> Vector{MatElem{C}}

Return the basis `basis(L)` of the Lie algebra `L` in the underlying matrix
representation.
"""
function matrix_repr_basis(L::LinearLieAlgebra{C}) where {C<:RingElement}
  return Vector{dense_matrix_type(C)}(L.basis)
end

@doc raw"""
    matrix_repr_basis(L::LinearLieAlgebra{C}, i::Int) -> MatElem{C}

Return the `i`-th element of the basis `basis(L)` of the Lie algebra `L` in the
underlying matrix representation.
"""
function matrix_repr_basis(L::LinearLieAlgebra{C}, i::Int) where {C<:RingElement}
  return (L.basis[i])::dense_matrix_type(C)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, ::MIME"text/plain", L::LinearLieAlgebra)
  io = pretty(io)
  println(io, _lie_algebra_type_to_string(get_attribute(L, :type, :unknown), L.n))
  println(io, Indent(), "of dimension $(dim(L))", Dedent())
  print(io, "over ")
  print(io, Lowercase(), coefficient_ring(L))
end

function Base.show(io::IO, L::LinearLieAlgebra)
  if get(io, :supercompact, false)
    print(io, _lie_algebra_type_to_compact_string(get_attribute(L, :type, :unknown), L.n))
  else
    io = pretty(io)
    print(
      io,
      _lie_algebra_type_to_string(get_attribute(L, :type, :unknown), L.n),
      " over ",
      Lowercase(),
    )
    print(IOContext(io, :supercompact => true), coefficient_ring(L))
  end
end

function _lie_algebra_type_to_string(type::Symbol, n::Int)
  if type == :general_linear
    return "General linear Lie algebra of degree $n"
  elseif type == :special_linear
    return "Special linear Lie algebra of degree $n"
  elseif type == :special_orthogonal
    return "Special orthogonal Lie algebra of degree $n"
  else
    return "Linear Lie algebra with $(n)x$(n) matrices"
  end
end

function _lie_algebra_type_to_compact_string(type::Symbol, n::Int)
  if type == :general_linear
    return "gl_$n"
  elseif type == :special_linear
    return "sl_$n"
  elseif type == :special_orthogonal
    return "so_$n"
  else
    return "Linear Lie algebra"
  end
end

function symbols(L::LinearLieAlgebra)
  return L.s
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

@doc raw"""
    coerce_to_lie_algebra_elem(L::LinearLieAlgebra{C}, x::MatElem{C}) -> LinearLieAlgebraElem{C}

Returns the element of `L` whose matrix representation corresponds to `x`.
If no such element exists, an error is thrown.
"""
function coerce_to_lie_algebra_elem(
  L::LinearLieAlgebra{C}, x::MatElem{C}
) where {C<:RingElement}
  @req size(x) == (L.n, L.n) "Invalid matrix dimensions."
  m = coefficient_vector(x, matrix_repr_basis(L))
  return L(m)
end

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

@doc raw"""
    matrix_repr(x::LinearLieAlgebraElem{C}) -> Mat{C}

Return the Lie algebra element `x` in the underlying matrix representation.
"""
function Generic.matrix_repr(x::LinearLieAlgebraElem)
  L = parent(x)
  return sum(
    c * b for (c, b) in zip(_matrix(x), matrix_repr_basis(L));
    init=zero_matrix(coefficient_ring(L), L.n, L.n),
  )
end

function bracket(
  x::LinearLieAlgebraElem{C}, y::LinearLieAlgebraElem{C}
) where {C<:RingElement}
  check_parent(x, y)
  L = parent(x)
  x_mat = matrix_repr(x)
  y_mat = matrix_repr(y)
  return coerce_to_lie_algebra_elem(L, x_mat * y_mat - y_mat * x_mat)
end

###############################################################################
#
#   Constructor
#
###############################################################################

@doc raw"""
    lie_algebra(R::Ring, n::Int, basis::Vector{<:MatElem{elem_type(R)}}, s::Vector{<:VarName}; cached::Bool) -> LinearLieAlgebra{elem_type(R)}

Construct the Lie algebra over the ring `R` with basis `basis` and basis element names
given by `s`. The basis elements must be square matrices of size `n`.
We require `basis` to be linearly independent, and to contain the Lie bracket of any
two basis elements in its span.

If `cached` is `true`, the constructed Lie algebra is cached.
"""
function lie_algebra(
  R::Ring,
  n::Int,
  basis::Vector{<:MatElem{C}},
  s::Vector{<:VarName};
  cached::Bool=true,
  check::Bool=true,
) where {C<:RingElement}
  return LinearLieAlgebra{elem_type(R)}(R, n, basis, Symbol.(s); cached)
end

function lie_algebra(
  basis::Vector{LinearLieAlgebraElem{C}}; check::Bool=true
) where {C<:RingElement}
  parent_L = parent(basis[1])
  @req all(parent(x) === parent_L for x in basis) "Elements not compatible."
  R = coefficient_ring(parent_L)
  n = parent_L.n
  s = map(AbstractAlgebra.obj_to_string, basis)
  return lie_algebra(R, n, matrix_repr.(basis), s; check)
end

@doc raw"""
    abelian_lie_algebra(R::Ring, n::Int) -> LinearLieAlgebra{elem_type(R)}
    abelian_lie_algebra(::Type{LinearLieAlgebra}, R::Ring, n::Int) -> LinearLieAlgebra{elem_type(R)}
    abelian_lie_algebra(::Type{AbstractLieAlgebra}, R::Ring, n::Int) -> AbstractLieAlgebra{elem_type(R)}

Return the abelian Lie algebra of dimension `n` over the ring `R`.
The first argument can be optionally provided to specify the type of the returned
Lie algebra.
"""
function abelian_lie_algebra(R::Ring, n::Int)
  @req n >= 0 "Dimension must be non-negative."
  return abelian_lie_algebra(LinearLieAlgebra, R, n)
end

function abelian_lie_algebra(::Type{T}, R::Ring, n::Int) where {T<:LinearLieAlgebra}
  @req n >= 0 "Dimension must be non-negative."
  basis = [(b = zero_matrix(R, n, n); b[i, i] = 1; b) for i in 1:n]
  s = ["x_$(i)" for i in 1:n]
  L = lie_algebra(R, n, basis, s; check=false)
  set_attribute!(L, :is_abelian => true)
  return L
end

@doc raw"""
    general_linear_lie_algebra(R::Ring, n::Int) -> LinearLieAlgebra{elem_type(R)}

Return the general linear Lie algebra $\mathfrak{gl}_n(R)$.

# Examples
```jldoctest
julia> L = general_linear_lie_algebra(QQ, 2)
General linear Lie algebra of degree 2
  of dimension 4
over rational field

julia> basis(L)
4-element Vector{LinearLieAlgebraElem{QQFieldElem}}:
 x_1_1
 x_1_2
 x_2_1
 x_2_2

julia> matrix_repr_basis(L)
4-element Vector{QQMatrix}:
 [1 0; 0 0]
 [0 1; 0 0]
 [0 0; 1 0]
 [0 0; 0 1]
```
"""
function general_linear_lie_algebra(R::Ring, n::Int)
  basis = [(b = zero_matrix(R, n, n); b[i, j] = 1; b) for i in 1:n for j in 1:n]
  s = ["x_$(i)_$(j)" for i in 1:n for j in 1:n]
  L = lie_algebra(R, n, basis, s; check=false)
  set_attribute!(L, :type => :general_linear)
  return L
end

@doc raw"""
    special_linear_lie_algebra(R::Ring, n::Int) -> LinearLieAlgebra{elem_type(R)}

Return the special linear Lie algebra $\mathfrak{sl}_n(R)$.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 2)
Special linear Lie algebra of degree 2
  of dimension 3
over rational field

julia> basis(L)
3-element Vector{LinearLieAlgebraElem{QQFieldElem}}:
 e_1_2
 f_1_2
 h_1

julia> matrix_repr_basis(L)
3-element Vector{QQMatrix}:
 [0 1; 0 0]
 [0 0; 1 0]
 [1 0; 0 -1]
```
"""
function special_linear_lie_algebra(R::Ring, n::Int)
  basis_e = [(b = zero_matrix(R, n, n); b[i, j] = 1; b) for i in 1:n for j in (i + 1):n]
  basis_f = [(b = zero_matrix(R, n, n); b[j, i] = 1; b) for i in 1:n for j in (i + 1):n]
  basis_h = [
    (b = zero_matrix(R, n, n); b[i, i] = 1; b[i + 1, i + 1] = -1; b) for i in 1:(n - 1)
  ]
  s_e = ["e_$(i)_$(j)" for i in 1:n for j in (i + 1):n]
  s_f = ["f_$(i)_$(j)" for i in 1:n for j in (i + 1):n]
  s_h = ["h_$(i)" for i in 1:(n - 1)]
  L = lie_algebra(R, n, [basis_e; basis_f; basis_h], [s_e; s_f; s_h]; check=false)
  set_attribute!(L, :type => :special_linear)
  return L
end

@doc raw"""
    special_orthogonal_lie_algebra(R::Ring, n::Int) -> LinearLieAlgebra{elem_type(R)}

Return the special orthogonal Lie algebra $\mathfrak{so}_n(R)$.

# Examples
```jldoctest
julia> L = special_orthogonal_lie_algebra(QQ, 3)
Special orthogonal Lie algebra of degree 3
  of dimension 3
over rational field

julia> basis(L)
3-element Vector{LinearLieAlgebraElem{QQFieldElem}}:
 x_1_2
 x_1_3
 x_2_3

julia> matrix_repr_basis(L)
3-element Vector{QQMatrix}:
 [0 1 0; -1 0 0; 0 0 0]
 [0 0 1; 0 0 0; -1 0 0]
 [0 0 0; 0 0 1; 0 -1 0]
```
"""
function special_orthogonal_lie_algebra(R::Ring, n::Int)
  basis = [
    (b = zero_matrix(R, n, n); b[i, j] = 1; b[j, i] = -1; b) for i in 1:n for j in (i + 1):n
  ]
  s = ["x_$(i)_$(j)" for i in 1:n for j in (i + 1):n]
  L = lie_algebra(R, n, basis, s; check=false)
  set_attribute!(L, :type => :special_orthogonal)
  return L
end
