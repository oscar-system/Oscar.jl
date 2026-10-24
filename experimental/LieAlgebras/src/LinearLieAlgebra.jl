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
  elseif type == :symplectic
    return "sp_$n"
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
    (L::LinearLieAlgebra{C})(x::MatElem{C}; check::Bool=true) -> LinearLieAlgebraElem{C}

Return the element of `L` whose matrix representation corresponds to `x`.

If no such element exists, an error is thrown.
With `check=false` this is not verified and `x` is assumed to lie in the span of `matrix_repr_basis(L)`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 2);

julia> x = L(matrix(QQ, [1 2; 3 -1]))
2*e_1_2 + 3*f_1_2 + h_1

julia> matrix_repr(x)
[1    2]
[3   -1]
```
"""
function (L::LinearLieAlgebra{C})(x::MatElem{C}; check::Bool=true) where {C<:FieldElem}
  m = _coefficient_row(L, x)
  check && @req _is_in_span(L, x, m) "Matrix is not in the span of the basis."
  return LieAlgebraElem(L, m)
end

@doc raw"""
    in(x::MatElem{C}, L::LinearLieAlgebra{C}) -> Bool

Return `true` if `x` is the matrix representation of an element of `L`, `false` otherwise.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 2);

julia> matrix(QQ, [1 2; 3 -1]) in L
true

julia> identity_matrix(QQ, 2) in L
false
```
"""
function Base.in(x::MatElem{C}, L::LinearLieAlgebra{C}) where {C<:FieldElem}
  return _is_in_span(L, x, _coefficient_row(L, x))
end

# Coefficient row of `x` with respect to `matrix_repr_basis(L)`, assuming `x` lies in their span.
function _coefficient_row(L::LinearLieAlgebra{C}, x::MatElem{C}) where {C<:FieldElem}
  @req base_ring(x) == coefficient_ring(L) "Incompatible base rings."
  @req size(x) == (L.n, L.n) "Invalid matrix dimensions."
  pivots, lgs_pivots_inv, _ = _matrix_repr_basis_pivot_data(L)
  v = matrix(
    coefficient_ring(L), 1, dim(L),
    C[x[div(k - 1, L.n) + 1, rem(k - 1, L.n) + 1] for k in pivots],
  )
  return v * lgs_pivots_inv
end

# Whether `x` equals the matrix with coefficient row `m = _coefficient_row(L, x)`.
function _is_in_span(
  L::LinearLieAlgebra{C}, x::MatElem{C}, m::MatElem{C}
) where {C<:FieldElem}
  _, _, lgs = _matrix_repr_basis_pivot_data(L)
  return m * lgs == _vec_row(x)
end

# Flatten `x` into a `1 x (nrows * ncols)` matrix.
function _vec_row(x::MatElem)
  return matrix(base_ring(x), 1, nrows(x) * ncols(x), _vec(transpose(x)))
end

# `lgs` is the `dim(L) x n^2` matrix with the flattened basis matrices as
# rows, `pivots` are `dim(L)` column indices such that `lgs[:, pivots]` is
# invertible, and `lgs_pivots_inv` is the inverse of that submatrix. Hence the
# coefficient vector of a matrix `x` in the span of the basis is
# `_vec_row(x)[:, pivots] * lgs_pivots_inv`.
@attr Tuple{Vector{Int},dense_matrix_type(C),dense_matrix_type(C)} function _matrix_repr_basis_pivot_data(
  L::LinearLieAlgebra{C}
) where {C<:FieldElem}
  lgs = zero_matrix(coefficient_ring(L), dim(L), L.n^2)
  for (k, b) in enumerate(matrix_repr_basis(L))
    lgs[k:k, :] = _vec_row(b)
  end
  r, ech = rref(lgs)
  @req r == dim(L) "Basis matrices are linearly dependent."
  pivots = Int[findfirst(j -> !is_zero_entry(ech, i, j), 1:ncols(ech)) for i in 1:r]
  return (pivots, inv(lgs[:, pivots]), lgs)
end

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

@doc raw"""
    matrix_repr(x::LinearLieAlgebraElem{C}; copy::Bool=true) -> MatElem{C}

Return the Lie algebra element `x` in the underlying matrix representation.

The matrix is cached in `x`.
With `copy=false` the cached matrix itself is returned and must not be mutated.
"""
function Generic.matrix_repr(
  x::LinearLieAlgebraElem{C}; copy::Bool=true
) where {C<:FieldElem}
  if !isdefined(x, :matrix_repr)
    L = parent(x)
    mat = zero_matrix(coefficient_ring(L), L.n, L.n)
    tmp = zero(mat)
    for (c, b) in zip(coefficients(x), matrix_repr_basis(L))
      iszero(c) && continue
      mat = addmul!(mat, b, c, tmp)
    end
    x.matrix_repr = mat
  end
  mat = (x.matrix_repr)::dense_matrix_type(C)
  return copy ? deepcopy(mat) : mat
end

function bracket(
  x::LinearLieAlgebraElem{C}, y::LinearLieAlgebraElem{C}
) where {C<:FieldElem}
  check_parent(x, y)
  L = parent(x)
  x_mat = matrix_repr(x; copy=false)
  y_mat = matrix_repr(y; copy=false)
  z_mat = x_mat * y_mat - y_mat * x_mat
  z = L(z_mat; check=false)
  z.matrix_repr = z_mat
  return z
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
  return lie_algebra(R, L.n, basis, symbols(L); is_basis=true, check=false)
end

###############################################################################
#
#   Constructor
#
###############################################################################

@doc raw"""
    lie_algebra(R::Field, n::Int, mats::Vector{<:MatElem}, s=:x; is_basis::Bool=false, check::Bool=true) -> LinearLieAlgebra{elem_type(R)}

Construct the Lie algebra generated by the $n \times n$ matrices `mats` over `R` under the
commutator bracket $[x, y] = xy - yx$, i.e. the smallest Lie subalgebra of
$\mathfrak{gl}_n(R)$ containing `mats`.

With `is_basis=true`, `mats` is used as basis instead, and `check=true` verifies that
`mats` is linearly independent and closed under the bracket.
Otherwise, `check` is currently ignored.

The basis elements are named by `s`, which is either a vector with one name per basis
element or a prefix such as `:x`, giving the names `x_1`, `x_2`, ...

# Examples
```jldoctest
julia> e = matrix(QQ, [0 0 1; 0 0 0; 0 0 0]);

julia> f = matrix(QQ, [0 0 0; 0 0 0; 1 0 0]);

julia> h = matrix(QQ, [1 0 0; 0 0 0; 0 0 -1]);

julia> L = lie_algebra(QQ, 3, [e, f, h], [:e, :f, :h]; is_basis=true)
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
  mats::Vector{<:MatElem},
  s::Union{VarName,Vector{<:VarName}}=:x;
  is_basis::Bool=false,
  check::Bool=true,
)
  @req all(g -> base_ring(g) === R, mats) "Matrices must be over the same field"
  @req all(g -> size(g) == (n, n), mats) "Matrices must be square and of the same size"

  if is_basis
    basis_mats = mats
  else
    gl = general_linear_lie_algebra(R, n)
    S = sub(gl, [gl(g; check=false) for g in mats])
    basis_mats = dense_matrix_type(R)[matrix_repr(b) for b in basis(S)]
    check = false # a basis of a subalgebra by construction
  end
  return LinearLieAlgebra{elem_type(R)}(
    R, n, basis_mats, _basis_symbols(s, length(basis_mats)); check
  )
end

_basis_symbols(s::VarName, k::Int) = [Symbol(s, :_, i) for i in 1:k]
_basis_symbols(s::Vector{<:VarName}, ::Int) = Symbol.(s)

function lie_algebra(
  elems::Vector{LinearLieAlgebraElem{C}}; is_basis::Bool=false, check::Bool=true
) where {C<:FieldElem}
  @req !isempty(elems) "List of elements must not be empty, or provide the Lie algebra as first argument"
  return lie_algebra(parent(elems[1]), elems; is_basis, check)
end

function lie_algebra(
  L::LinearLieAlgebra{C},
  elems::Vector{LinearLieAlgebraElem{C}};
  is_basis::Bool=false,
  check::Bool=true,
) where {C<:FieldElem}
  @req all(parent(x) === L for x in elems) "Elements not compatible."
  R = coefficient_ring(L)
  if is_basis
    basis_elems = elems
  else
    basis_elems = basis(sub(L, elems))
    check = false # a basis of a subalgebra by construction
  end
  s = map(AbstractAlgebra.obj_to_string_wrt_times, basis_elems)
  return lie_algebra(R, L.n, matrix_repr.(basis_elems), s; is_basis=true, check)
end

@doc raw"""
    lie_algebra(gens::Vector{<:MatElem{C}}, s=:x; is_basis::Bool=false, check::Bool=true) -> LinearLieAlgebra{C}

Construct the Lie algebra generated by the square matrices `gens` under the
commutator bracket $[x, y] = xy - yx$, i.e. the smallest Lie subalgebra of
$\mathfrak{gl}_n(R)$ containing `gens`.

All matrices in `gens` must be of the same size $n \times n$ over the same field $R$.
See [`lie_algebra(::Field, ::Int, ::Vector{<:MatElem})`](@ref) for `s`, `is_basis` and
`check`.

# Examples
```jldoctest
julia> e = matrix(QQ, [0 1 0; 0 0 1; 0 0 0]);

julia> f = matrix(QQ, [0 0 0; 2 0 0; 0 2 0]);

julia> L = lie_algebra([e, f])
Linear Lie algebra with 3x3 matrices
  of dimension 3
over rational field

julia> matrix_repr_basis(L)
3-element Vector{QQMatrix}:
 [1 0 0; 0 0 0; 0 0 -1]
 [0 1 0; 0 0 1; 0 0 0]
 [0 0 0; 1 0 0; 0 1 0]
```
"""
function lie_algebra(
  gens::Vector{<:MatElem{C}},
  s::Union{VarName,Vector{<:VarName}}=:x;
  is_basis::Bool=false,
  check::Bool=true,
) where {C<:FieldElem}
  @req !isempty(gens) "List of matrices must not be empty, or provide the field and size as first arguments"
  return lie_algebra(base_ring(gens[1]), nrows(gens[1]), gens, s; is_basis, check)
end

@doc raw"""
    lie_algebra(R::Field, gens::Vector{<:MatElem}, s=:x; is_basis::Bool=false, check::Bool=true) -> LinearLieAlgebra{elem_type(R)}

Construct the Lie algebra over `R` generated by the square matrices `gens`, whose entries are
mapped into `R`.

See [`lie_algebra(::Field, ::Int, ::Vector{<:MatElem})`](@ref) for `s`, `is_basis` and
`check`.

# Examples
```jldoctest
julia> e = matrix(ZZ, [0 1 0; 0 0 1; 0 0 0]);

julia> f = matrix(ZZ, [0 0 0; 2 0 0; 0 2 0]);

julia> dim(lie_algebra(QQ, [e, f]))
3

julia> dim(lie_algebra(GF(2), [e, f]))
1
```
"""
function lie_algebra(
  R::Field,
  gens::Vector{<:MatElem},
  s::Union{VarName,Vector{<:VarName}}=:x;
  is_basis::Bool=false,
  check::Bool=true,
)
  @req !isempty(gens) "List of matrices must not be empty"
  if all(g -> base_ring(g) === R, gens)
    return lie_algebra(R, nrows(gens[1]), gens, s; is_basis, check)
  end
  mapped_gens = [map_entries(R, g) for g in gens]
  return lie_algebra(R, nrows(mapped_gens[1]), mapped_gens, s; is_basis, check)
end

function abelian_lie_algebra(::Type{T}, R::Field, n::Int) where {T<:LinearLieAlgebra}
  @req n >= 0 "Dimension must be non-negative."
  basis = [(b = zero_matrix(R, n, n); b[i, i] = 1; b) for i in 1:n]
  L = lie_algebra(R, n, basis; is_basis=true, check=false)
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
  L = lie_algebra(R, n, basis, s; is_basis=true, check=false)
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
  L = lie_algebra(
    R, n, [basis_e; basis_f; basis_h], [s_e; s_f; s_h]; is_basis=true, check=false
  )
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

If `gram` is not provided, for $n = 2k$ the form defined by $\begin{pmatrix} 0 & I_k \\ I_k & 0 \end{pmatrix}$
is used, and for $n = 2k + 1$ the form defined by $\begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & I_k \\ 0 & I_k & 0 \end{pmatrix}$.

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
  L = lie_algebra(R, n, basis; is_basis=true, check=false)
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

If `gram` is not provided, for $n = 2k$ the form defined by $\begin{pmatrix} 0 & I_k \\ -I_k & 0 \end{pmatrix}$
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
  L = lie_algebra(R, n, basis; is_basis=true, check=false)
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
