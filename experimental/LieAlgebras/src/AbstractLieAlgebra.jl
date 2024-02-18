@attributes mutable struct AbstractLieAlgebra{C<:FieldElem} <: LieAlgebra{C}
  R::Field
  dim::Int
  struct_consts::Matrix{SRow{C}}
  s::Vector{Symbol}

  # only set if known
  root_system::RootSystem

  function AbstractLieAlgebra{C}(
    R::Field,
    struct_consts::Matrix{SRow{C}},
    s::Vector{Symbol};
    cached::Bool=true,
    check::Bool=true,
  ) where {C<:FieldElem}
    return get_cached!(
      AbstractLieAlgebraDict, (R, struct_consts, s), cached
    ) do
      (n1, n2) = size(struct_consts)
      @req n1 == n2 "Invalid structure constants dimensions."
      dimL = n1
      @req length(s) == dimL "Invalid number of basis element names."
      if check
        @req all(
          r -> all(e -> parent(last(e)) === R, r), struct_consts
        ) "Invalid structure constants."
        @req all(
          iszero, struct_consts[i, i][k] for i in 1:dimL, k in 1:dimL
        ) "Not anti-symmetric."
        @req all(
          iszero,
          struct_consts[i, j][k] + struct_consts[j, i][k] for i in 1:dimL, j in 1:dimL,
          k in 1:dimL
        ) "Not anti-symmetric."
        @req all(
          iszero,
          sum(
            struct_consts[i, j][k] * struct_consts[k, l][m] +
            struct_consts[j, l][k] * struct_consts[k, i][m] +
            struct_consts[l, i][k] * struct_consts[k, j][m] for k in 1:dimL
          ) for i in 1:dimL, j in 1:dimL, l in 1:dimL, m in 1:dimL
        ) "Jacobi identity does not hold."
      end
      new{C}(R, dimL, struct_consts, s)
    end::AbstractLieAlgebra{C}
  end
end

const AbstractLieAlgebraDict = CacheDictType{
  Tuple{Field,Matrix{SRow},Vector{Symbol}},AbstractLieAlgebra
}()

struct AbstractLieAlgebraElem{C<:FieldElem} <: LieAlgebraElem{C}
  parent::AbstractLieAlgebra{C}
  mat::MatElem{C}
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{AbstractLieAlgebraElem{C}}) where {C<:FieldElem} = AbstractLieAlgebra{C}

elem_type(::Type{AbstractLieAlgebra{C}}) where {C<:FieldElem} = AbstractLieAlgebraElem{C}

parent(x::AbstractLieAlgebraElem) = x.parent

coefficient_ring(L::AbstractLieAlgebra{C}) where {C<:FieldElem} = L.R::parent_type(C)

dim(L::AbstractLieAlgebra) = L.dim

###############################################################################
#
#   Root system getters
#
###############################################################################

has_root_system(L::LieAlgebra) = isdefined(L, :root_system)

function root_system(L::LieAlgebra)
  @req has_root_system(L) "No root system known."
  return L.root_system
end

@doc raw"""
    chevalley_basis(L::AbstractLieAlgebra{C}) -> NTuple{3,Vector{AbstractLieAlgebraElem{C}}}

Return the Chevalley basis of the Lie algebra `L` in three vectors, stating first the positive root vectors, 
then the negative root vectors, and finally the basis of the Cartan subalgebra. The order of root vectors corresponds
to the order of the roots in the root system.
"""
function chevalley_basis(L::AbstractLieAlgebra)
  @req has_root_system(L) "No root system known."

  npos = n_positive_roots(root_system(L))
  b = basis(L)
  # root vectors
  r_plus = b[1:npos]
  r_minus = b[(npos + 1):(2 * npos)]
  # basis for cartan algebra
  h = b[(2 * npos + 1):dim(L)]
  return (r_plus, r_minus, h)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, ::MIME"text/plain", L::AbstractLieAlgebra)
  io = pretty(io)
  println(io, "Abstract Lie algebra")
  println(io, Indent(), "of dimension $(dim(L))", Dedent())
  print(io, "over ")
  print(io, Lowercase(), coefficient_ring(L))
end

function Base.show(io::IO, L::AbstractLieAlgebra)
  if get(io, :supercompact, false)
    print(io, "Abstract Lie algebra")
  else
    io = pretty(io)
    print(io, "Abstract Lie algebra over ", Lowercase())
    print(IOContext(io, :supercompact => true), coefficient_ring(L))
  end
end

function symbols(L::AbstractLieAlgebra)
  return L.s
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

# no special ones

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function bracket(
  x::AbstractLieAlgebraElem{C}, y::AbstractLieAlgebraElem{C}
) where {C<:FieldElem}
  check_parent(x, y)
  L = parent(x)
  mat = sum(
    cxi * cyj * L.struct_consts[i, j] for (i, cxi) in enumerate(coefficients(x)),
    (j, cyj) in enumerate(coefficients(y));
    init=sparse_row(coefficient_ring(L)),
  )
  return L(mat)
end

###############################################################################
#
#   Properties
#
###############################################################################

function is_abelian(L::AbstractLieAlgebra)
  return all(e -> iszero(length(e)), L.struct_consts)
end

###############################################################################
#
#   Constructor
#
###############################################################################

@doc raw"""
    lie_algebra(R::Field, struct_consts::Matrix{SRow{elem_type(R)}}, s::Vector{<:VarName}; cached::Bool, check::Bool) -> AbstractLieAlgebra{elem_type(R)}

Construct the Lie algebra over the ring `R` with structure constants `struct_consts`
and with basis element names `s`.

The Lie bracket on the newly constructed Lie algebra `L` is determined by the structure
constants in `struct_consts` as follows: let $x_i$ denote the $i$-th standard basis vector
of `L`. Then the entry `struct_consts[i,j][k]` is a scalar $a_{i,j,k}$
such that $[x_i, x_j] = \sum_k a_{i,j,k} x_k$.

* `s`: A vector of basis element names. This is 
  `[Symbol("x_$i") for i in 1:size(struct_consts, 1)]` by default.
* `cached`: If `true`, cache the result. This is `true` by default.
* `check`: If `true`, check that the structure constants are anti-symmetric and
  satisfy the Jacobi identity. This is `true` by default.
"""
function lie_algebra(
  R::Field,
  struct_consts::Matrix{SRow{C}},
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:size(struct_consts, 1)];
  cached::Bool=true,
  check::Bool=true,
) where {C<:FieldElem}
  @req C == elem_type(R) "Invalid coefficient type."
  return AbstractLieAlgebra{elem_type(R)}(R, struct_consts, Symbol.(s); cached, check)
end

@doc raw"""
    lie_algebra(R::Field, struct_consts::Array{elem_type(R),3}, s::Vector{<:VarName}; cached::Bool, check::Bool) -> AbstractLieAlgebra{elem_type(R)}

Construct the Lie algebra over the ring `R` with structure constants `struct_consts`
and with basis element names `s`.

The Lie bracket on the newly constructed Lie algebra `L` is determined by the structure
constants in `struct_consts` as follows: let $x_i$ denote the $i$-th standard basis vector
of `L`. Then the entry `struct_consts[i,j,k]` is a scalar $a_{i,j,k}$
such that $[x_i, x_j] = \sum_k a_{i,j,k} x_k$.

* `s`: A vector of basis element names. This is
  `[Symbol("x_$i") for i in 1:size(struct_consts, 1)]` by default.
* `cached`: If `true`, cache the result. This is `true` by default.
* `check`: If `true`, check that the structure constants are anti-symmetric and
  satisfy the Jacobi identity. This is `true` by default.

# Examples
```jldoctest
julia> struct_consts = zeros(QQ, 3, 3, 3);

julia> struct_consts[1, 2, 3] = QQ(1);

julia> struct_consts[2, 1, 3] = QQ(-1);

julia> struct_consts[3, 1, 1] = QQ(2);

julia> struct_consts[1, 3, 1] = QQ(-2);

julia> struct_consts[3, 2, 2] = QQ(-2);

julia> struct_consts[2, 3, 2] = QQ(2);

julia> sl2 = lie_algebra(QQ, struct_consts, ["e", "f", "h"])
Abstract Lie algebra
  of dimension 3
over rational field

julia> e, f, h = basis(sl2)
3-element Vector{AbstractLieAlgebraElem{QQFieldElem}}:
 e
 f
 h

julia> e * f
h

julia> h * e
2*e

julia> h * f
-2*f
```
"""
function lie_algebra(
  R::Field,
  struct_consts::Array{C,3},
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:size(struct_consts, 1)];
  cached::Bool=true,
  check::Bool=true,
) where {C<:FieldElem}
  @req C == elem_type(R) "Invalid coefficient type."
  struct_consts2 = Matrix{SRow{elem_type(R)}}(
    undef, size(struct_consts, 1), size(struct_consts, 2)
  )
  for i in axes(struct_consts, 1), j in axes(struct_consts, 2)
    struct_consts2[i, j] = sparse_row(
      R, collect(axes(struct_consts, 3)), struct_consts[i, j, :]
    )
  end

  return AbstractLieAlgebra{elem_type(R)}(R, struct_consts2, Symbol.(s); cached, check)
end

function lie_algebra(
  basis::Vector{AbstractLieAlgebraElem{C}}; check::Bool=true
) where {C<:FieldElem}
  parent_L = parent(basis[1])
  @req all(parent(x) === parent_L for x in basis) "Elements not compatible."
  R = coefficient_ring(parent_L)
  basis_matrix = if length(basis) == 0
    matrix(R, 0, dim(L), C[])
  else
    matrix(R, [coefficients(b) for b in basis])
  end
  struct_consts = Matrix{SRow{elem_type(R)}}(undef, length(basis), length(basis))
  for (i, bi) in enumerate(basis), (j, bj) in enumerate(basis)
    fl, row = can_solve_with_solution(basis_matrix, _matrix(bi * bj); side=:left)
    @req fl "Not closed under the bracket."
    struct_consts[i, j] = sparse_row(row)
  end
  s = map(AbstractAlgebra.obj_to_string, basis)
  return lie_algebra(R, struct_consts, s; check)
end

@doc raw"""
    lie_algebra(R::Field, dynkin::Tuple{Char,Int}; cached::Bool) -> AbstractLieAlgebra{elem_type(R)}

Construct the simple Lie algebra over the ring `R` with Dynkin type given by `dynkin`.
The internally used basis of this Lie algebra is the Chevalley basis.

If `cached` is `true`, the constructed Lie algebra is cached.
"""
function lie_algebra(R::Field, S::Symbol, n::Int; cached::Bool=true)
  rs = root_system(S, n)
  cm = cartan_matrix(rs)
  @req is_cartan_matrix(cm; generalized=false) "The type does not correspond to a classical root system"

  npos = n_positive_roots(rs)
  nsimp = n_simple_roots(rs)
  n = 2 * npos + nsimp

  #=
  struct_consts = Matrix{SRow{elem_type(R)}}(undef, n, n)
  for i in 1:npos, j in 1:npos
    # [x_i, x_j]
    fl, k = is_positive_root_with_index(positive_root(rs, i) + positive_root(rs, j))
    struct_consts[i, j] = fl ? sparse_row(R, [k], [1]) : sparse_row(R)
    # [x_i, y_j] = δ_ij h_i
    struct_consts[i, npos + j] = i == j ? sparse_row(R, [2 * npos + i], [1]) : sparse_row(R)
    # [y_j, x_i] = -[x_i, y_j]
    struct_consts[npos + j, i] = -struct_consts[i, npos + j]
    # [y_i, y_j]
    fl, k = is_negative_root_with_index(negative_root(rs, i) + negative_root(rs, j))
    struct_consts[npos + i, npos + j] = fl ? sparse_row(R, [npos + k], [1]) : sparse_row(R)
  end
  for i in 1:nsimp, j in 1:npos
    # [h_i, x_j] = <α_j, α_i> x_j
    struct_consts[2 * npos + i, j] = sparse_row(R, [j], [cm[j, i]])
    # [h_i, y_j] = - <α_j, α_i> y_j
    struct_consts[2 * npos + i, npos + j] = sparse_row(R, [npos + j], [-cm[j, i]])
    # [x_j, h_i] = -[h_i, x_j]
    struct_consts[j, 2 * npos + i] = -struct_consts[2 * npos + i, j]
    # [y_j, h_i] = -[h_i, y_j]
    struct_consts[npos + j, 2 * npos + i] = -struct_consts[2 * npos + i, npos + j]
  end
  for i in 1:nsimp, j in 1:nsimp
    # [h_i, h_j] = 0
    struct_consts[2 * npos + i, 2 * npos + j] = sparse_row(R)
  end

  s = [
    [Symbol("x_$i") for i in 1:npos]
    [Symbol("y_$i") for i in 1:npos]
    [Symbol("h_$i") for i in 1:nsimp]
  ]

  L = lie_algebra(R, struct_consts, s; cached, check=true) # TODO: remove check
  =#

  # start temporary workaround # TODO: reenable code above
  type = only(root_system_type(rs))
  if type == (:F, 4)
    # GAP uses a non-canonical order of simple roots for F4. Until we have our own implementation, we disable it to avoid confusion.
    error("Not implemented for F4.")
  end
  coeffs_iso = inv(Oscar.iso_oscar_gap(R))
  LG = GAP.Globals.SimpleLieAlgebra(GAP.Obj(string(type[1])), type[2], domain(coeffs_iso))
  @req GAPWrap.Dimension(LG) == n "Dimension mismatch. Something went wrong."
  s = [
    [Symbol("x_$i") for i in 1:npos]
    [Symbol("y_$i") for i in 1:npos]
    [Symbol("h_$i") for i in 1:nsimp]
  ]
  L = codomain(
    _iso_gap_oscar_abstract_lie_algebra(LG, s; coeffs_iso, cached)
  )::AbstractLieAlgebra{elem_type(R)}
  # end temporary workaround

  set_attribute!(L, :is_simple, true)
  return L
end

function abelian_lie_algebra(::Type{T}, R::Field, n::Int) where {T<:AbstractLieAlgebra}
  @req n >= 0 "Dimension must be non-negative."
  basis = [(b = zero_matrix(R, n, n); b[i, i] = 1; b) for i in 1:n]
  s = ["x_$(i)" for i in 1:n]
  L = lie_algebra(R, n, basis, s; check=false)

  struct_consts = Matrix{SRow{elem_type(R)}}(undef, n, n)
  for i in axes(struct_consts, 1), j in axes(struct_consts, 2)
    struct_consts[i, j] = sparse_row(R)
  end

  L = lie_algebra(R, struct_consts, s)
  set_attribute!(L, :is_abelian => true)
  return L
end
