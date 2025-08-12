###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LinearLieAlgebraElem{C}}) where {C<:FieldElem} = LinearLieAlgebra{C}

elem_type(::Type{LinearLieAlgebra{C}}) where {C<:FieldElem} = LinearLieAlgebraElem{C}

parent(x::LinearLieAlgebraElem) = x.parent

coefficient_ring(L::LinearLieAlgebra{C}) where {C<:FieldElem} = L.R::parent_type(C)

vector_space_dim(L::LinearLieAlgebra) = L.dim

@doc raw"""
    matrix_repr_basis(L::LinearLieAlgebra{C}) -> Vector{MatElem{C}}

Return the basis `basis(L)` of the Lie algebra `L` in the underlying matrix
representation.
"""
function matrix_repr_basis(L::LinearLieAlgebra{C}) where {C<:FieldElem}
  return L.basis::Vector{dense_matrix_type(C)}
end

@doc raw"""
    matrix_repr_basis(L::LinearLieAlgebra{C}, i::Int) -> MatElem{C}

Return the `i`-th element of the basis `basis(L)` of the Lie algebra `L` in the
underlying matrix representation.
"""
function matrix_repr_basis(L::LinearLieAlgebra{C}, i::Int) where {C<:FieldElem}
  return matrix_repr_basis(L)[i]
end

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, mime::MIME"text/plain", L::LinearLieAlgebra)
  @show_name(io, L)
  @show_special(io, mime, L)
  io = pretty(io)
  type_string = _lie_algebra_type_to_string(get_attribute(L, :type, :unknown), L.n)
  if !isnothing(type_string)
    println(io, type_string)
  else
    println(io, "Linear Lie algebra with $(L.n)x$(L.n) matrices")
    if has_root_system(L)
      rs = root_system(L)
      if has_root_system_type(rs)
        type, ord = root_system_type_with_ordering(rs)
        print(io, Indent(), "of type ", _root_system_type_string(type))
        if !issorted(ord)
          print(io, " (non-canonical ordering)")
        end
        println(io, Dedent())
      end
    end
  end
  println(io, Indent(), "of dimension ", dim(L), Dedent())
  print(io, "over ", Lowercase(), coefficient_ring(L))
end

function Base.show(io::IO, L::LinearLieAlgebra)
  @show_name(io, L)
  @show_special(io, L)
  if is_terse(io)
    type_string_compact = _lie_algebra_type_to_compact_string(
      get_attribute(L, :type, :unknown), L.n
    )
    if !isnothing(type_string_compact)
      print(io, type_string_compact)
    else
      print(io, "Linear Lie algebra")
    end
  else
    io = pretty(io)
    type_string = _lie_algebra_type_to_string(get_attribute(L, :type, :unknown), L.n)
    if !isnothing(type_string)
      print(io, type_string)
    else
      print(io, "Linear Lie algebra with $(L.n)x$(L.n) matrices")
      if has_root_system(L)
        rs = root_system(L)
        if has_root_system_type(rs)
          type, ord = root_system_type_with_ordering(rs)
          print(io, " of type ", _root_system_type_string(type))
          if !issorted(ord)
            print(io, " (non-canonical ordering)")
          end
        end
      end
    end
    print(terse(io), " over ", Lowercase(), coefficient_ring(L))
  end
end

function _lie_algebra_type_to_string(type::Symbol, n::Int)
  if type == :general_linear
    return "General linear Lie algebra of degree $n"
  elseif type == :special_linear
    return "Special linear Lie algebra of degree $n"
  elseif type == :special_orthogonal
    return "Special orthogonal Lie algebra of degree $n"
  elseif type == :symplectic
    return "Symplectic Lie algebra of degree $n"
  end
  return nothing
end

function _lie_algebra_type_to_compact_string(type::Symbol, n::Int)
  if type == :general_linear
    return "gl_$n"
  elseif type == :special_linear
    return "sl_$n"
  elseif type == :special_orthogonal
    return "so_$n"
  end
  return nothing
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

Return the element of `L` whose matrix representation corresponds to `x`.
If no such element exists, an error is thrown.
"""
function coerce_to_lie_algebra_elem(
  L::LinearLieAlgebra{C}, x::MatElem{C}
) where {C<:FieldElem}
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
  mat = zero_matrix(coefficient_ring(L), L.n, L.n)
  tmp = zero(mat)
  for (c, b) in zip(coefficients(x), matrix_repr_basis(L))
    mat = addmul!(mat, b, c, tmp)
  end
  return mat
end

function bracket(
  x::LinearLieAlgebraElem{C}, y::LinearLieAlgebraElem{C}
) where {C<:FieldElem}
  check_parent(x, y)
  L = parent(x)
  x_mat = matrix_repr(x)
  y_mat = matrix_repr(y)
  return coerce_to_lie_algebra_elem(L, x_mat * y_mat - y_mat * x_mat)
end

###############################################################################
#
#   Root system getters
#
###############################################################################

has_root_system(L::LinearLieAlgebra) = isdefined(L, :root_system)

function root_system(L::LinearLieAlgebra)
  assure_root_system(L)
  return L.root_system
end

function chevalley_basis(L::LinearLieAlgebra)
  assure_root_system(L)
  return L.chevalley_basis::NTuple{3,Vector{elem_type(L)}}
end

function set_root_system_and_chevalley_basis!(
  L::LinearLieAlgebra{C}, R::RootSystem, chev::NTuple{3,Vector{LinearLieAlgebraElem{C}}}
) where {C<:FieldElem}
  L.root_system = R
  L.chevalley_basis = chev
end

###############################################################################
#
#   change_base_ring
#
###############################################################################

function change_base_ring(R::Field, L::LinearLieAlgebra)
  basis = map(b -> change_base_ring(R, b), matrix_repr_basis(L))
  return lie_algebra(R, L.n, basis, symbols(L); check=false)
end

###############################################################################
#
#   Constructor
#
###############################################################################

@doc raw"""
    lie_algebra(R::Field, n::Int, basis::Vector{<:MatElem{elem_type(R)}}, s::Vector{<:VarName}; check::Bool=true) -> LinearLieAlgebra{elem_type(R)}

Construct the Lie algebra over the field `R` with basis `basis` and basis element names
given by `s`. The basis elements must be square matrices of size `n`.

We require `basis` to be linearly independent, and to contain the Lie bracket of any
two basis elements in its span (this is currently not checked).
Setting `check=false` disables these checks (once they are in place).

# Examples
```jldoctest
julia> e = matrix(QQ, [0 0 1; 0 0 0; 0 0 0]);

julia> f = matrix(QQ, [0 0 0; 0 0 0; 1 0 0]);

julia> h = matrix(QQ, [1 0 0; 0 0 0; 0 0 -1]);

julia> L = lie_algebra(QQ, 3, [e, f, h], [:e, :f, :h])
Linear Lie algebra with 3x3 matrices
  of dimension 3
over rational field

julia> root_system(L);

julia> L
Linear Lie algebra with 3x3 matrices
  of type A1
  of dimension 3
over rational field

julia> basis(L)
3-element Vector{LinearLieAlgebraElem{QQFieldElem}}:
 e
 f
 h

julia> matrix_repr_basis(L)
3-element Vector{QQMatrix}:
 [0 0 1; 0 0 0; 0 0 0]
 [0 0 0; 0 0 0; 1 0 0]
 [1 0 0; 0 0 0; 0 0 -1]
```
"""
function lie_algebra(
  R::Field,
  n::Int,
  basis::Vector{<:MatElem{C}},
  s::Vector{<:VarName};
  check::Bool=true,
) where {C<:FieldElem}
  return LinearLieAlgebra{elem_type(R)}(R, n, basis, Symbol.(s); check)
end

function lie_algebra(
  basis::Vector{LinearLieAlgebraElem{C}}; check::Bool=true
) where {C<:FieldElem}
  @req !isempty(basis) "Basis must not be empty, or provide the Lie algebra as first argument"
  return lie_algebra(parent(basis[1]), basis; check)
end

function lie_algebra(
  L::LinearLieAlgebra{C}, basis::Vector{LinearLieAlgebraElem{C}}; check::Bool=true
) where {C<:FieldElem}
  @req all(parent(x) === L for x in basis) "Elements not compatible."
  R = coefficient_ring(L)
  s = map(AbstractAlgebra.obj_to_string_wrt_times, basis)
  return lie_algebra(R, L.n, matrix_repr.(basis), s; check)
end

function abelian_lie_algebra(::Type{T}, R::Field, n::Int) where {T<:LinearLieAlgebra}
  @req n >= 0 "Dimension must be non-negative."
  basis = [(b = zero_matrix(R, n, n); b[i, i] = 1; b) for i in 1:n]
  s = ["x_$(i)" for i in 1:n]
  L = lie_algebra(R, n, basis, s; check=false)
  set_attribute!(L, :is_abelian => true)
  return L
end

@doc raw"""
    general_linear_lie_algebra(R::Field, n::Int) -> LinearLieAlgebra{elem_type(R)}

Return the general linear Lie algebra $\mathfrak{gl}_n(R)$,
i.e., the Lie algebra of all $n \times n$ matrices over the field `R`.

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
function general_linear_lie_algebra(R::Field, n::Int)
  basis = [(b = zero_matrix(R, n, n); b[i, j] = 1; b) for i in 1:n for j in 1:n]
  s = ["x_$(i)_$(j)" for i in 1:n for j in 1:n]
  L = lie_algebra(R, n, basis, s; check=false)
  set_attribute!(L, :type => :general_linear)
  return L
end

@doc raw"""
    special_linear_lie_algebra(R::Field, n::Int) -> LinearLieAlgebra{elem_type(R)}

Return the special linear Lie algebra $\mathfrak{sl}_n(R)$,
i.e., the Lie algebra of all $n \times n$ matrices over the field `R` with trace zero.

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
function special_linear_lie_algebra(R::Field, n::Int)
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

function _lie_algebra_basis_from_form(R::Field, n::Int, form::MatElem)
  invform = inv(form)
  eqs = zero_matrix(R, n^2, n^2)
  for i in 1:n, j in 1:n
    x = zero_matrix(R, n, n)
    x[i, j] = 1
    eqs[(i - 1) * n + j, :] = _vec(x + invform * transpose(x) * form)
  end
  ker = kernel(eqs)
  rref!(ker) # we cannot assume anything about the kernel, but want to have a consistent output
  dim = nrows(ker)
  basis = [zero_matrix(R, n, n) for _ in 1:dim]
  for i in 1:n
    for k in 1:dim
      basis[k][i, 1:n] = ker[k, (i - 1) * n .+ (1:n)]
    end
  end
  return basis
end

@doc raw"""
    special_orthogonal_lie_algebra(R::Field, n::Int) -> LinearLieAlgebra{elem_type(R)}
    special_orthogonal_lie_algebra(R::Field, n::Int, gram::MatElem) -> LinearLieAlgebra{elem_type(R)}
    special_orthogonal_lie_algebra(R::Field, n::Int, gram::Matrix) -> LinearLieAlgebra{elem_type(R)}

Return the special orthogonal Lie algebra $\mathfrak{so}_n(R)$.

Given a non-degenerate symmetric bilinear form $f$ via its Gram matrix `gram`,
$\mathfrak{so}_n(R)$ is the Lie algebra of all $n \times n$ matrices $x$ over the field `R`
such that $f(xv, w) = -f(v, xw)$ for all $v, w \in R^n$.

If `gram` is not provided, for $n = 2k$ the form defined by $\begin{matrix} 0 & I_k \\ -I_k & 0 \end{matrix}$
is used, and for $n = 2k + 1$ the form defined by $\begin{matrix} 1 & 0 & 0 \\ 0 & 0 I_k \\ 0 & I_k & 0 \end{matrix}$.

# Examples
```jldoctest
julia> L1 = special_orthogonal_lie_algebra(QQ, 4)
Special orthogonal Lie algebra of degree 4
  of dimension 6
over rational field

julia> basis(L1)
6-element Vector{LinearLieAlgebraElem{QQFieldElem}}:
 x_1
 x_2
 x_3
 x_4
 x_5
 x_6

julia> matrix_repr_basis(L1)
6-element Vector{QQMatrix}:
 [1 0 0 0; 0 0 0 0; 0 0 -1 0; 0 0 0 0]
 [0 1 0 0; 0 0 0 0; 0 0 0 0; 0 0 -1 0]
 [0 0 0 1; 0 0 -1 0; 0 0 0 0; 0 0 0 0]
 [0 0 0 0; 1 0 0 0; 0 0 0 -1; 0 0 0 0]
 [0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 -1]
 [0 0 0 0; 0 0 0 0; 0 1 0 0; -1 0 0 0]

julia> L2 = special_orthogonal_lie_algebra(QQ, 3, identity_matrix(QQ, 3))
Special orthogonal Lie algebra of degree 3
  of dimension 3
over rational field

julia> basis(L2)
3-element Vector{LinearLieAlgebraElem{QQFieldElem}}:
 x_1
 x_2
 x_3

julia> matrix_repr_basis(L2)
3-element Vector{QQMatrix}:
 [0 1 0; -1 0 0; 0 0 0]
 [0 0 1; 0 0 0; -1 0 0]
 [0 0 0; 0 0 1; 0 -1 0]
```
"""
special_orthogonal_lie_algebra

function special_orthogonal_lie_algebra(R::Field, n::Int, gram::MatElem)
  form = map_entries(R, gram)
  @req size(form) == (n, n) "Invalid matrix dimensions"
  @req is_symmetric(form) "Bilinear form must be symmetric"
  @req is_invertible(form) "Bilinear form must be non-degenerate"
  basis = _lie_algebra_basis_from_form(R, n, form)
  dim = length(basis)
  @assert characteristic(R) != 0 || dim == div(n^2 - n, 2)
  s = ["x_$(i)" for i in 1:dim]
  L = lie_algebra(R, n, basis, s; check=false)
  set_attribute!(L, :type => :special_orthogonal, :form => form)
  return L
end

function special_orthogonal_lie_algebra(R::Field, n::Int, gram::Matrix)
  return special_orthogonal_lie_algebra(R, n, matrix(R, gram))
end

function special_orthogonal_lie_algebra(R::Field, n::Int)
  if is_even(n)
    k = div(n, 2)
    form = zero_matrix(R, n, n)
    form[1:k, k .+ (1:k)] = identity_matrix(R, k)
    form[k .+ (1:k), 1:k] = identity_matrix(R, k)
  else
    k = div(n - 1, 2)
    form = zero_matrix(R, n, n)
    form[1, 1] = one(R)
    form[1 .+ (1:k), 1 + k .+ (1:k)] = identity_matrix(R, k)
    form[1 + k .+ (1:k), 1 .+ (1:k)] = identity_matrix(R, k)
  end
  return special_orthogonal_lie_algebra(R, n, form)
end

@doc raw"""
    symplectic_lie_algebra(R::Field, n::Int) -> LinearLieAlgebra{elem_type(R)}
    symplectic_lie_algebra(R::Field, n::Int, gram::MatElem) -> LinearLieAlgebra{elem_type(R)}
    symplectic_lie_algebra(R::Field, n::Int, gram::Matrix) -> LinearLieAlgebra{elem_type(R)}

Return the symplectic Lie algebra $\mathfrak{sp}_n(R)$.

Given a non-degenerate skew-symmetric bilinear form $f$ via its Gram matrix `gram`,
$\mathfrak{sp}_n(R)$ is the Lie algebra of all $n \times n$ matrices $x$ over the field `R`
such that $f(xv, w) = -f(v, xw)$ for all $v, w \in R^n$.

If `gram` is not provided, for $n = 2k$ the form defined by $\begin{matrix} 0 & I_k \\ -I_k & 0 \end{matrix}$
is used.
For odd $n$ there is no non-degenerate skew-symmetric bilinear form on $R^n$.

# Examples
```jldoctest
julia> L = symplectic_lie_algebra(QQ, 4)
Symplectic Lie algebra of degree 4
  of dimension 10
over rational field

julia> basis(L)
10-element Vector{LinearLieAlgebraElem{QQFieldElem}}:
 x_1
 x_2
 x_3
 x_4
 x_5
 x_6
 x_7
 x_8
 x_9
 x_10

julia> matrix_repr_basis(L)
10-element Vector{QQMatrix}:
 [1 0 0 0; 0 0 0 0; 0 0 -1 0; 0 0 0 0]
 [0 1 0 0; 0 0 0 0; 0 0 0 0; 0 0 -1 0]
 [0 0 1 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
 [0 0 0 1; 0 0 1 0; 0 0 0 0; 0 0 0 0]
 [0 0 0 0; 1 0 0 0; 0 0 0 -1; 0 0 0 0]
 [0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 -1]
 [0 0 0 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
 [0 0 0 0; 0 0 0 0; 1 0 0 0; 0 0 0 0]
 [0 0 0 0; 0 0 0 0; 0 1 0 0; 1 0 0 0]
 [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 1 0 0]
```
"""
symplectic_lie_algebra

function symplectic_lie_algebra(R::Field, n::Int, gram::MatElem)
  form = map_entries(R, gram)
  @req size(form) == (n, n) "Invalid matrix dimensions"
  @req is_skew_symmetric(form) "Bilinear form must be skew-symmetric"
  @req is_even(n) && is_invertible(form) "Bilinear form must be non-degenerate"
  basis = _lie_algebra_basis_from_form(R, n, form)
  dim = length(basis)
  @assert characteristic(R) != 0 || dim == div(n^2 + n, 2)
  s = ["x_$(i)" for i in 1:dim]
  L = lie_algebra(R, n, basis, s; check=false)
  set_attribute!(L, :type => :symplectic, :form => form)
  return L
end

function symplectic_lie_algebra(R::Field, n::Int, gram::Matrix)
  return symplectic_lie_algebra(R, n, matrix(R, gram))
end

function symplectic_lie_algebra(R::Field, n::Int)
  @req is_even(n) "Dimension must be even"
  k = div(n, 2)
  form = zero_matrix(R, n, n)
  form[1:k, k .+ (1:k)] = identity_matrix(R, k)
  form[k .+ (1:k), 1:k] = -identity_matrix(R, k)
  return symplectic_lie_algebra(R, n, form)
end
